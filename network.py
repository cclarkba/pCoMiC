####### NEURAL NETWORK #######

import torch.nn as nn
import torch.nn.functional as F
import torch
import numpy as np
import matplotlib.pyplot as plt
from data_loader import MutDataset, MutDataset_Unsuper

##### TODO Look into multi-head attention network, I think it is similar to the multi-stream thing you're already doing
#####   but could be some minute differences, transformer is likely a good route to take since the features we are learning
#####   are so wildly different


class Transformer(nn.Module):
    # For variable input lengths, could be interesting to use this in tandem with other properties and get rid of first two networks below
    # Alternatively, if this works well, just use only this and makes life a lot easier
    # Haven't debugged this at all so will probably need to look over after finished processing training data
    # network = Transformer(1480, 128, 8, 256, 2)

    def __init__(self, input_length, embedding_size, num_heads, ff_dim, num_transformer_layers):
        super(Transformer, self).__init__()

        self.input_length = input_length

        self.embedding = nn.Linear(20, embedding_size)


        self.transformer_layers = nn.ModuleList([
            nn.TransformerEncoderLayer(
                d_model=embedding_size,
                nhead=num_heads,
                dim_feedforward=ff_dim,
                activation='relu'
            )
            for _ in range(num_transformer_layers)
        ])

        self.global_avg_pool = nn.AdaptiveAvgPool1d(1)
        self.fc1 = nn.Linear(embedding_size, 30)
        self.fc2 = nn.Linear(30,20)
        self.output_layer = nn.Linear(20,5)
        self.sigmoid = nn.Sigmoid()
    
    def forward(self, x):
        x = x.float()
        x = self.embedding(x)
        for layer in self.transformer_layers:
            x = layer(x)
        x = x.permute(0, 2, 1)
        x = self.global_avg_pool(x).squeeze(-1)
        x = self.fc1(x)
        x = self.fc2(x)
        x = self.output_layer(x)
        x = self.sigmoid(x)
        return x
    
def semisup_train(network, train_dl, test_dl, unlab_dl, optimizer, criterion, num_epochs, alpha=0.9):
    network.train()
    
    losses = []
    train_accs = []
    test_accs = []

    def calc_acc(dataloader):
        correct = 0
        total = 0

        with torch.no_grad():
            for data,target in dataloader:
                data = data.float()
                outputs = network(data)
                predicted = (outputs >= 0.5).float()
                correct += (predicted == target).all(dim=1).sum().item()
                total += len(target)
        accuracy = (correct/total)*100
        return accuracy
    
    for epoch in range(num_epochs):
        total_loss = 0.0

        for (lab_data, lab_target), unlab_data in zip(train_dl, unlab_dl):
            # Forward pass through labeled
            lab_data = lab_data.float()
            lab_target = lab_target.float()
            lab_out = network(lab_data)
            lab_loss = criterion(lab_out, lab_target)

            # Foward through unlab
            unlab_data = unlab_data.float()
            unlab_out = network(unlab_data)
            pseudo_lab = (unlab_out>0.5).float()
            consistency_loss = F.mse_loss(unlab_out, pseudo_lab)

            combined_loss = alpha*lab_loss + (1-alpha)*consistency_loss

            optimizer.zero_grad()
            combined_loss.backward()
            optimizer.step()

            total_loss += combined_loss.item()
        
        total_loss = total_loss / len(train_dl)
        train_accs.append(calc_acc(train_dl))
        test_accs.append(calc_acc(test_dl))
        losses.append(total_loss)

        print(f'Epoch {epoch + 1}/{num_epochs}: Loss = {total_loss:.4f}')
        print(f'Train acc: {round(train_accs[-1],2)}%, Test acc: {round(test_accs[-1],2)}%')
    return train_accs, test_accs, losses

def train(network, dataloader, dataloader_test, optimizer, criterion, num_epochs=10):
    network.train()

    losses = []
    train_accs = []
    test_accs = []
    def calc_acc(dataloader):
        correct = 0
        total = 0

        with torch.no_grad():
            for data,target in dataloader:
                data = data.float()
                outputs = network(data)
                predicted = (outputs >= 0.5).float()
                correct += (predicted == target).all(dim=1).sum().item()
                total += len(target)
        accuracy = (correct/total)*100
        return accuracy

    for epoch in range(num_epochs):
        total_loss = 0.0
        correct = 0
        total = 0
        for data,target in dataloader:
            data = data.float()
            optimizer.zero_grad()
            outputs = network(data)
            target = target.float()
            loss = criterion(outputs, target)
            loss.backward()
            optimizer.step()
            total_loss += loss.item()

            predicted = (outputs >= 0.5).float()
            correct += (predicted == target).all(dim=1).sum().item()
            total += len(target)

        avg_loss = total_loss / len(dataloader)
        train_acc = (correct/total) *100
        test_acc = calc_acc(dataloader_test)
        losses.append(avg_loss)
        train_accs.append(train_acc)
        test_accs.append(test_acc)
        print(f'Epoch {epoch+1} Loss: {avg_loss:.4f}')
        print(f'Train acc: {train_acc:.4f}%, Test_acc: {test_acc:.4f}%')
        network.train()
    return train_accs, test_accs, losses

# Loading in all the data
def split_labeled_data(file, test_size):
    from sklearn.model_selection import train_test_split
    data = torch.load(file)
    dataset = data.dataset
    targets = data.labels

    # Need to remove singleton classes (ie classes only represented once) and ensure they are placed in test set
    # Count occurrences of each class using a dictionary
    class_counts = {}
    for idx, target in enumerate(targets):
        class_item = tuple(target.detach().numpy())
        if class_counts.get(class_item):
            class_counts[class_item].append(idx)
        else:
            class_counts[class_item] = [idx]

    # Identify, remove and store singleton classes
    singleton_indices = []
    for class_label, indices in class_counts.items():
        if len(indices) == 1:  # Singleton class
            singleton_indices.extend(indices)

    # Filter out the singleton samples
    remaining_indices = [i for i in range(len(dataset)) if i not in singleton_indices]
    remaining_targets = targets[remaining_indices]

    # Create a subset of the remaining data
    dataset = torch.utils.data.Subset(dataset, remaining_indices)
    targets = remaining_targets

    # Step 3: Proceed with train-test split
    test_size = 0.2  # 20% for testing, 80% for training
    random_state = 42  # For reproducibility

    # Perform train-test split
    train_indices, test_indices = train_test_split(
        range(len(dataset)), test_size=test_size, stratify=targets, random_state=random_state
    )
    test_indices.extend(singleton_indices)
    # Create training and testing subsets
    train_subset = torch.utils.data.Subset(data, train_indices)
    test_subset = torch.utils.data.Subset(data, test_indices)

    return train_subset, test_subset

def create_data_loader(train_set, test_set, unlabel_file=None):
    label_dl_train = torch.utils.data.DataLoader(train_set, batch_size=len(train_set)//4, shuffle=True)
    label_dl_test = torch.utils.data.DataLoader(test_set, batch_size=len(test_set)//4, shuffle=False)
    if unlabel_file:
        unlabeled = torch.load(unlabel_file)
        unlabel_dl = torch.utils.data.DataLoader(unlabeled, batch_size=len(unlabeled)//4, shuffle=True)
        return label_dl_train, label_dl_test, unlabel_dl
    else:
        return label_dl_train, label_dl_test

def test_accuracy(network, label_dl_train, label_dl_test):
    correct = 0
    total_len = 0
    network.eval()
    threshold = 0.5
    for data, target in label_dl_train:
        out = network(data)
        bin_pred = (out > threshold).float()
        match_mask = (bin_pred == target).float()
        matches = match_mask*bin_pred
        per_sample = torch.sum(matches, dim=1)
        correct += torch.sum(per_sample).item()
        total_len += len(data)
    train_accuracy = correct / total_len

    correct = 0
    total_len = 0
    for data, target in label_dl_test:
        out = network(data)
        bin_pred = (out > threshold).float()
        match_mask = (bin_pred == target).float()
        matches = match_mask*bin_pred
        per_sample = torch.sum(matches, dim=1)
        correct += torch.sum(per_sample).item()
        total_len += len(data)
    test_accuracy = correct / total_len
    return train_accuracy, test_accuracy

"""#### Setting parameters for the network itself ####
input_length = 1480  # Set the input length (sequence length) according to your data
embedding_size = 64  # Size of the embedding vectors
num_heads = 4  # Number of attention heads
ff_dim = 128  # Dimensionality of the feed-forward network
num_transformer_layers = 2  # Number of transformer encoder layers

#### Creating network and related parameters ####
network = Transformer(input_length, embedding_size, num_heads, ff_dim, num_transformer_layers)
criterion = nn.BCEWithLogitsLoss() # maybe try just BCELoss
# criterion = nn.BCELoss()
optimizer = torch.optim.Adam(network.parameters(), lr=0.001)

#### Loading and processing data ####
train_subset, test_subset = split_labeled_data("C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/dataset_split.pt", test_size=0.2)
label_dl_train, label_dl_test, unlabel_dl = create_data_loader(train_subset, test_subset, 
                                                               "C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/dataset_unsuper_trimmed.pt")

#### Two different training methods ####
# train_accs, test_accs, losses = train(network, label_dl_train, label_dl_test, optimizer, criterion, num_epochs=100)
train_accs, test_accs, label_losses, total_losses = semisup_train(network, label_dl_train, label_dl_test, unlabel_dl, optimizer, criterion, num_epochs=100)

#### Saving trained network ####
# torch.save(network, "C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/transformer_superv2.pth")
torch.save(network, "C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/transformer_semisuperv2.pth")

#### Load trained network ####
# network = torch.load("C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/transformer_semisuper.pth")
# network = torch.load("C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/transformer_superv2.pth")

#### Test Accuracy ####
# test_accuracy(network, label_dl_train, label_dl_test)
# import matplotlib.pyplot as plt

# fig, ax = plt.subplots(1,2)
# ax[0].plot(train_accs, 'r-')
# ax[0].plot(test_accs, 'b-')
# ax[0].set_title('Accuracies')

# # lab_loss_copy = label_losses.numpy()
# # tot_loss_copy = total_losses.numpy()

# ax[1].plot(losses,'r-')
# # ax[1].plot(tot_loss_copy, 'b-')
# ax[1].set_title('Losses')

# plt.show()"""


##### Going to make a new fully connected network which takes in ratio of normal/mutated length, hydropathy, charge and hbond scores
##### Hopefully this will at least classify some correctly

class Network(nn.Module):
    def __init__(self):
        super(Network, self).__init__()
        self.fc1 = nn.Linear(4, 50)
        self.fc2 = nn.Linear(50,50)
        self.fc3 = nn.Linear(50,50)
        self.output = nn.Linear(50,5)
        self.relu = nn.ReLU()
        self.sigmoid = nn.Sigmoid()

    
    def forward(self, x):
        x = self.fc1(x)
        x = self.relu(x)
        x = self.fc2(x)
        x = self.relu(x)
        x = self.fc3(x)
        x = self.relu(x)
        x = self.output(x)
        x = self.sigmoid(x)
        return x



# network = Transformer(1480, 128, 8, 256, 1, 2)
# network = Network()
# criterion = nn.BCELoss() # maybe try just BCELoss
# optimizer = torch.optim.NAdam(network.parameters(), lr=0.001)

# #### Loading and processing data ####
train_subset, test_subset = split_labeled_data("C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/dataset_char_scores_super.pt", test_size=0.2)
label_dl_train, label_dl_test, unlab_dl = create_data_loader(train_subset, test_subset, 
                                                   "C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/dataset_char_scores_unsuper.pt")

# train_subset, test_subset = split_labeled_data("C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/dataset_split.pt", test_size=0.2)
# label_dl_train, label_dl_test, unlab_dl = create_data_loader(train_subset, test_subset, 
#                                                    "C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/dataset_unsuper_trimmed.pt")


# # #### Two different training methods ####
# # # train_accs, test_accs, losses = train(network, label_dl_train, label_dl_test, optimizer, criterion, num_epochs=200)
# train_accs, test_accs, losses = semisup_train(network, label_dl_train, label_dl_test, unlab_dl, optimizer, criterion, num_epochs=200)

# fig,ax = plt.subplots(1,2, figsize=(11,4))
# ax[0].plot(train_accs,'r-', label='Training')
# ax[0].plot(test_accs,'b-', label='Testing')
# ax[0].set_xlabel('Epochs')
# ax[0].set_ylabel('Accuracy (%)')
# ax[0].set_title('Accuracy over Epochs for Dense Network')

# ax[1].plot(losses)
# ax[1].set_xlabel('Epochs')
# ax[1].set_ylabel('BCELoss')
# ax[1].set_title('Loss over Epochs for Dense Network')
# plt.show()

### TODO Make plots showing accuracies of each class for each network (4 plots, make sure to show benign)
###     additional note: correct classification should be deemed if any of the expected classes are represented
        # ie, if target is [0 1 1 0 0], [0 1 0 0 0] and [0 1 1 0 0] or any other form also including a 1 in the second or third position is correct


# train_dense(network, label_dl_train, label_dl_test, criterion, optimizer, num_epochs=200)

# alpha = np.arange(.0,1.0,0.1)
# final_train = [0.0, 0.0, 0.0, 0.4634, 0.8293, 0.8049, 0.8049, 0.8537, 0.8049, 0.8293]
# final_test = [0.0, 0.0, 0.0, 0.2308, 0.3846, 0.3077, 0.3077, 0.5385, 0.4615, 0.5385]

# i,j,k = train_dense_semi(network, label_dl_train, label_dl_test, unlab_dl, criterion, optimizer, alpha=0.5, num_epochs=200)

# plt.plot(alpha, final_train, 'r-', label='Train')
# plt.plot(alpha, final_test, 'b-', label='Test')
# plt.xlabel('alpha')
# plt.ylabel('Accuracy (%)')
# plt.title('Effects of Varying Alpha on Accuracy')
# plt.legend()
# plt.show()


#### Saving trained network ####
# torch.save(network, "C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/TEST_transformer_semisuper.pth")
# torch.save(network, "C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/FINAL_transformer_super.pth")
# torch.save(network, "C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/FINAL_dense_semisuper.pth")
# torch.save(network, "C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/FINAL_dense_super.pth")

#### Load trained network ####
# network = torch.load("C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/FINAL_transformer_semisuper.pth")
# network = torch.load("C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/FINAL_transformer_super.pth")
network = torch.load("C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/FINAL_dense_semisuper.pth")
# network = torch.load("C:/Users/musta/OneDrive/Desktop/pCoMiC/test_data/FINAL_dense_super.pth")

classes = ['Production', 'Processing', 'Gating', 'Conducting', 'Insufficient', 'Benign']

probs_train = {
    'Production': [0,0],
    'Processing': [0,0],
    'Gating': [0,0],
    'Conducting': [0,0],
    'Insufficient': [0,0],
    'Benign': [0,0]
}

threshold = 0.6
network.eval()
for data, target in label_dl_train:
    data = data.float()
    out = network(data)
    for d,t in zip(out, target):
        d = d.detach().numpy()
        t = t.detach().numpy()

        highest = np.argmax(d)
        targ_indices = [i for i,x in enumerate(t) if x == 1]

        update_both = False
        if highest in targ_indices:
            update_both = True
        
        for idx in targ_indices:
            class_add = classes[idx]
            probs_train[class_add][1] += 1
            if update_both:
                probs_train[class_add][0] += 1
        if not targ_indices:
            probs_train['Benign'][1] += 1
            if np.all([x for x in d if x < 0.6]):
                probs_train['Benign'][0] += 1    

probs_test = {
    'Production': [0,0],
    'Processing': [0,0],
    'Gating': [0,0],
    'Conducting': [0,0],
    'Insufficient': [0,0],
    'Benign': [0,0]
}

threshold = 0.6
network.eval()
for data, target in label_dl_test:
    data = data.float()
    out = network(data)
    for d,t in zip(out, target):
        d = d.detach().numpy()
        t = t.detach().numpy()

        highest = np.argmax(d)
        targ_indices = [i for i,x in enumerate(t) if x == 1]

        update_both = False
        if highest in targ_indices:
            update_both = True
        
        for idx in targ_indices:
            class_add = classes[idx]
            probs_test[class_add][1] += 1
            if update_both:
                probs_test[class_add][0] += 1
        if not targ_indices:
            probs_test['Benign'][1] += 1
            if np.all([x for x in d if x < 0.6]):
                probs_test['Benign'][0] += 1

cls = list(probs_train.keys())
train_probs = []
for p in probs_train.values():
    if p[1] != 0:
        train_probs.append((p[0]/p[1]))
    else:
        train_probs.append(0)

test_probs = []
for p in probs_test.values():
    if p[1] != 0:
        test_probs.append((p[0]/p[1]))
    else:
        test_probs.append(0)

print(probs_train)
print(probs_test)
fig, ax = plt.subplots()
bar_container = ax.bar(cls, test_probs)
ax.set(ylabel='Correct Identification (%)', title='Correct Identification among Classes for Dense Net (Semi-Sup, Test)', ylim=(0, 1.1))
ax.bar_label(bar_container, fmt='{:,.2f}')
plt.show()

