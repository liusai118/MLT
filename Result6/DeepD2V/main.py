import torch
import torch.nn as nn
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, StratifiedKFold
from torch.utils.data import DataLoader
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, confusion_matrix, roc_curve, auc
import matplotlib.pyplot as plt
import seaborn as sns
import json
import os

# 假设 MyDataSet 和 BCL_Network 已经定义在 myDataSet 和 models 模块中
import myDataSet as ms
import load_data as ld
from models import BCL_Network


def create_folder(path, fold):
    fold_path = os.path.join(path, f"fold_{fold}")
    os.makedirs(fold_path, exist_ok=True)
    return fold_path


def path():
    return "D:/Programming/python/PycharmProjects/ProteinDNABinding/"


# ROC曲线绘制函数
def draw_ROC_Curve(probabilities, y_true, save_path):
    # 使用预测概率计算 FPR 和 TPR
    fpr, tpr, _ = roc_curve(y_true, probabilities)
    roc_auc = auc(fpr, tpr)

    # 绘制 ROC 曲线
    plt.figure()
    plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (area = {roc_auc:.2f})')
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC) Curve')
    plt.legend(loc="lower right")

    # 保存图像为 SVG 格式
    plt.savefig(save_path, format='svg')
    plt.close()
    return roc_auc


# 训练函数，返回最优模型的路径
def train(myDataLoader, validate_DataLoader, path, fold):
    best_f1 = 0
    train_metrics = []
    validate_metrics = []

    fold_path = create_folder(path, fold)

    for epoch in range(Epoch):
        for step, (x, y) in enumerate(myDataLoader):
            model.train()
            x, y = x.to(device), y.to(device)
            output = model(x)
            loss = loss_func(output, y)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        # 计算训练集指标
        train_f1, _, train_metrics_values = validate(myDataLoader, epoch, fold_path, is_training=True)
        train_metrics.append(train_metrics_values)

        # 计算验证集指标
        validate_f1, test_loss, validate_metrics_values = validate(validate_DataLoader, epoch, fold_path)
        validate_metrics.append(validate_metrics_values)

        if validate_f1 > best_f1:
            best_f1 = validate_f1
            best_model_name = os.path.join(fold_path, f'validate_params_epoch_{epoch}.pkl')
            torch.save(model.state_dict(), best_model_name)

    # 保存训练和验证集的指标为 JSON 文件
    with open(os.path.join(fold_path, 'train_metrics.json'), 'w') as f:
        json.dump(train_metrics, f, indent=4)
    with open(os.path.join(fold_path, 'validate_metrics.json'), 'w') as f:
        json.dump(validate_metrics, f, indent=4)

    scheduler.step(test_loss)
    print(f"Best model saved at: {best_model_name}")
    return best_model_name


# 验证函数
def validate(myDataLoader, epoch, fold_path, is_training=False):
    output_list = []
    correct_list = []
    test_loss = 0
    for step, (x, y) in enumerate(myDataLoader):
        model.eval()
        x, y = x.to(device), y.to(device)
        output = model(x)
        loss = loss_func(output, y)
        test_loss += float(loss)
        output_list += output.cpu().detach().numpy().tolist()
        correct_list += y.cpu().detach().numpy().tolist()

    y_pred = (np.array(output_list) > 0.5).astype(int).flatten()
    y_true = np.array(correct_list).flatten()

    accuracy = accuracy_score(y_true, y_pred)
    precision = precision_score(y_true, y_pred, average="weighted")
    recall = recall_score(y_true, y_pred, average="weighted")
    F1 = f1_score(y_true, y_pred, average="weighted")
    test_loss /= len(myDataLoader)

    # 保存混淆矩阵图，改为 SVG 格式
    cm = confusion_matrix(y_true, y_pred)
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt="d", cmap="Blues")
    plt.xlabel("Predicted Label")
    plt.ylabel("True Label")
    plt.title(f"Confusion Matrix - {'Training' if is_training else 'Validation'} Epoch {epoch}")
    plt.savefig(os.path.join(fold_path, f'confusion_matrix_epoch_{epoch}.svg'), format='svg')  # 改为 SVG
    plt.close()

    # 将指标保存为字典
    metrics = {
        "epoch": epoch,
        "f1": F1,
        "precision": precision,
        "accuracy": accuracy,
        "recall": recall
    }

    return F1, test_loss, metrics


# 评估函数
def evaluate(myDataLoader, path, fold, model_path):
    fold_path = create_folder(path, fold)
    model.load_state_dict(torch.load(model_path, map_location=device))
    model.eval()

    probabilities = []
    correct_list = []
    for step, (x, y) in enumerate(myDataLoader):
        x, y = x.to(device), y.to(device)
        with torch.no_grad():
            output = model(x)
            probabilities += output.cpu().numpy().tolist()  # 确保使用概率
        correct_list += y.cpu().numpy().tolist()

    y_true = np.array(correct_list).flatten()
    probabilities = np.array(probabilities).flatten()

    # 计算评价指标
    y_pred = (probabilities > 0.5).astype(int)
    accuracy = accuracy_score(y_true, y_pred)
    precision = precision_score(y_true, y_pred, average="weighted")
    recall = recall_score(y_true, y_pred, average="weighted")
    F1 = f1_score(y_true, y_pred, average="weighted")

    # 保存 ROC 曲线并计算 AUC
    roc_auc = draw_ROC_Curve(probabilities, y_true, os.path.join(fold_path, f'evaluation_roc_curve.svg'))

    return roc_auc, precision, F1


# 数据加载函数
def getDataSet(train_index, validate_index):
    x_train = X.iloc[train_index]
    y_train = y.iloc[train_index]
    x_validate = X.iloc[validate_index]
    y_validate = y.iloc[validate_index]
    x_train_, x_validate_, y_train_, y_validate_ = train_test_split(
        x_train, y_train, test_size=0.125, stratify=y_train, random_state=1)
    x_train_ = x_train_.reset_index(drop=True)
    x_validate_ = x_validate_.reset_index(drop=True)
    x_validate = x_validate.reset_index(drop=True)
    y_train_ = y_train_.reset_index(drop=True)
    y_validate_ = y_validate_.reset_index(drop=True)
    y_validate = y_validate.reset_index(drop=True)

    train_DataSet = ms.MyDataSet(input=x_train_, label=y_train_)
    validate_DataSet = ms.MyDataSet(input=x_validate_, label=y_validate_)
    validate_Full_DataSet = ms.MyDataSet(input=x_validate, label=y_validate)
    train_DataLoader = DataLoader(dataset=train_DataSet, batch_size=Batch_Size, shuffle=True)
    validate_DataLoader = DataLoader(dataset=validate_DataSet, batch_size=test_Batch_Size, shuffle=True)
    validate_Full_DataLoader = DataLoader(dataset=validate_Full_DataSet, batch_size=test_Batch_Size, shuffle=False)
    return train_DataLoader, validate_DataLoader, validate_Full_DataLoader


# 主程序
if __name__ == '__main__':
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    Batch_Size = 64
    test_Batch_Size = 256
    LR = 3e-4
    Epoch = 5
    K_Fold = 5
    file_list = ld.create_list(ld.dataDir)
    file_list.sort()
    file_list = ['NRF22']
    num_gpus = torch.cuda.device_count()
    device_ids = list(range(num_gpus))

    for path in file_list:
        all_data = pd.read_csv(ld.dataDir + path + '/all_data.txt', sep='\t')
        X = all_data.iloc[:, 0]  # 假设第一列是 sequence
        y = all_data.iloc[:, 1]  # 假设第二列是 label
        kf = StratifiedKFold(n_splits=K_Fold, shuffle=True, random_state=1)
        fold = 1
        roc_total, pr_total, F1_total = 0, 0, 0

        for train_index, validate_index in kf.split(X, y):
            train_DataLoader, validate_DataLoader, validate_Full_DataLoader = getDataSet(train_index, validate_index)
            model = BCL_Network().to(device)
            if num_gpus > 1:
                model = nn.DataParallel(model, device_ids=device_ids)
            optimizer = torch.optim.Adam(model.parameters(), lr=LR)
            scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', factor=0.5, patience=5)
            loss_func = nn.BCELoss()

            best_model_name = train(train_DataLoader, validate_DataLoader, path, fold)
            roc_auc, precision, F1 = evaluate(validate_Full_DataLoader, path, fold, best_model_name)

            roc_total += roc_auc
            pr_total += precision
            F1_total += F1

            fold += 1
        roc_average = roc_total / K_Fold
        pr_average = pr_total / K_Fold
        f1_average = F1_total / K_Fold
        print(f"Path: {path}, Avg ROC: {roc_average}, Avg PR: {pr_average}, Avg F1: {f1_average}")
