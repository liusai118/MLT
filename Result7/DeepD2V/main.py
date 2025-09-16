import torch
import torch.nn as nn
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from torch.utils.data import Dataset, DataLoader
from sklearn.metrics import roc_auc_score, roc_curve, auc, accuracy_score, recall_score, f1_score, \
    precision_recall_curve
from sklearn.model_selection import StratifiedKFold
import matplotlib.pyplot as plt
from torch.utils.tensorboard import SummaryWriter
import os
# 假设 MyDataSet 和 BCL_Network 已经定义在 myDataSet 和 models 模块中
import myDataSet as ms
import load_data as ld
from models import BCL_Network
import util

from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

def path():
    return "/root/autodl-tmp/pn/"


# 训练函数，返回最优模型的路径
def train(myDataLoader, validate_DataLoader, path, fold, writer):
    best_roc = 0
    global_step = 0  # 用于记录步数
    train_f1_scores = []
    validate_f1_scores = []
    loss_values = []  # 用于记录每个step的loss值

    for epoch in range(Epoch):
        for step, (x, y) in enumerate(myDataLoader):
            model.train()
            x, y = x.to(device), y.to(device)
            output = model(x)
            loss = loss_func(output, y)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            # 记录每一步的loss值
            loss_values.append(loss.item())
            global_step += 1

        _, _, train_f1, _, _ = validate(myDataLoader, epoch, save_csv=False)
        train_f1_scores.append(train_f1)

        ROC, PR, validate_f1, test_loss, accuracy = validate(validate_DataLoader, epoch)
        validate_f1_scores.append(validate_f1)

        if ROC > best_roc:
            best_roc = ROC
            torch.save(model.state_dict(),
                       ld.modelDir + path + '/validate_params_' + str(fold) + '_' + str(epoch) + '.pkl')
            best_model_name = ld.modelDir + path + '/validate_params_' + str(fold) + '_' + str(epoch) + '.pkl'

    scheduler.step(test_loss)
    print(f"Best model saved at: {best_model_name}")

    # Save F1 scores for plotting
    f1_scores = {
        "train": train_f1_scores,
        "validate": validate_f1_scores
    }

    return best_model_name, f1_scores, loss_values


def validate(myDataLoader, epoch, save_csv=False, csv_path="validation_predictions.csv"):
    output_list = []
    output_result_list = []
    correct_list = []
    test_loss = 0
    for step, (x, y) in enumerate(myDataLoader):
        model.eval()
        x, y = x.to(device), y.to(device)
        output = model(x)
        loss = loss_func(output, y)
        test_loss += float(loss)
        output_list += output.cpu().detach().numpy().tolist()
        output = (output > 0.5).int()
        output_result_list += output.cpu().detach().numpy().tolist()
        correct_list += y.cpu().detach().numpy().tolist()

    y_pred = np.array(output_result_list).flatten()
    y_true = np.array(correct_list).flatten()
    accuracy = accuracy_score(y_true, y_pred)
    test_loss /= len(myDataLoader)
    print(f'Validation: Avg. loss: {test_loss:.4f}, Accuracy: {accuracy:.3f}')
    ROC, PR, F1 = util.draw_ROC_Curve(output_list, output_result_list, correct_list, path + '/' + 'test')

    if save_csv:
        predictions_df = pd.DataFrame({"True Label": y_true, "Predicted Label": y_pred})
        predictions_df.to_csv(csv_path, index=False)
        print(f"Validation predictions saved to {csv_path}")

    return ROC, PR, F1, test_loss, accuracy


def predict(model_path, test_loader):
    # 加载模型
    model = BCL_Network().to(device)
    model.load_state_dict(torch.load(model_path, map_location=device))
    model.eval()  # 切换到评估模式

    all_outputs = []
    all_labels = []

    with torch.no_grad():
        for inputs, labels in test_loader:
            inputs, labels = inputs.to(device), labels.to(device)
            outputs = model(inputs)
            # 将回归值转换为二进制标签
            binary_outputs = (outputs > 0.5).int()
            all_outputs.append(binary_outputs.cpu().numpy())
            all_labels.append(labels.cpu().numpy())

    # 确保预测结果和标签是1维数组
    return np.concatenate(all_outputs).flatten(), np.concatenate(all_labels).flatten()


from sklearn.metrics import (
    accuracy_score, f1_score, precision_score, recall_score,
    roc_curve, auc, confusion_matrix, ConfusionMatrixDisplay
)
import matplotlib.pyplot as plt
import numpy as np
import os

def evaluate(myDataLoader, path, fold, best_model_name):
    model.load_state_dict(torch.load(best_model_name, map_location=device))
    model.eval()

    output_list = []
    output_result_list = []
    correct_list = []

    for step, (x, y) in enumerate(myDataLoader):
        x, y = x.to(device), y.to(device)
        with torch.no_grad():
            output = model(x)
        output_list += output.cpu().numpy().tolist()
        binary_output = (output > 0.5).int()
        output_result_list += binary_output.cpu().numpy().tolist()
        correct_list += y.cpu().numpy().tolist()

    y_probs = np.array(output_list).flatten()
    y_pred = np.array(output_result_list).flatten()
    y_true = np.array(correct_list).flatten()

    # === 指标计算 ===
    accuracy = accuracy_score(y_true, y_pred)
    f1 = f1_score(y_true, y_pred)
    precision = precision_score(y_true, y_pred)
    recall = recall_score(y_true, y_pred)

    # === 图像保存路径 ===
    figure_dir = os.path.join("figure", path)
    os.makedirs(figure_dir, exist_ok=True)

    # === ROC 曲线 ===
    fpr, tpr, _ = roc_curve(y_true, y_probs)
    roc_auc = auc(fpr, tpr)

    plt.figure()
    plt.plot(fpr, tpr, label=f'ROC Curve (AUC = {roc_auc:.2f})')
    plt.plot([0, 1], [0, 1], 'k--', lw=1)
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'ROC Curve - Fold {fold}')
    plt.legend(loc='lower right')
    plt.grid(True)
    plt.savefig(f"{figure_dir}/roc_curve_fold_{fold}.svg")
    plt.close()

    # === 混淆矩阵 ===
    cm = confusion_matrix(y_true, y_pred)
    disp = ConfusionMatrixDisplay(confusion_matrix=cm)
    disp.plot(cmap='Blues', values_format='d')
    plt.title(f'Confusion Matrix - Fold {fold}')
    plt.savefig(f"{figure_dir}/confusion_matrix_fold_{fold}.svg")
    plt.close()

    # === 调用 util 中自定义曲线绘制（保留）===
    ROC, PR, _ = util.draw_ROC_Curve(output_list, output_result_list, correct_list,
                                     path + '/validate_params_' + str(fold))

    # === 返回所有必要结果 ===
    return ROC, PR, f1, accuracy, precision, recall, y_true, y_probs


def plot_roc_all_folds(fold_results, path):
    plt.figure(figsize=(8, 6))
    for result in fold_results:
        y_true = np.array(result["y_true"]).flatten()
        y_prob = np.array(result["y_prob"]).flatten()
        fold = result["fold"]

        fpr, tpr, _ = roc_curve(y_true, y_prob)
        roc_auc = auc(fpr, tpr)
        plt.plot(fpr, tpr, lw=2, label=f'Fold {fold} (AUC = {roc_auc:.2f})')

    plt.plot([0, 1], [0, 1], 'k--', lw=1)
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'ROC Curves Across All Folds')
    plt.legend(loc='lower right')
    plt.grid(True)

    out_dir = os.path.join("figure", path)
    os.makedirs(out_dir, exist_ok=True)
    plt.savefig(f"{out_dir}/roc_curve_all_folds.svg")
    plt.close()



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


def plot_loss_curve(loss_values, fold, smoothing_window=50):
    steps = range(1, len(loss_values) + 1)

    # 平滑损失值
    smoothed_loss_values = pd.Series(loss_values).rolling(window=smoothing_window).mean()

    plt.figure(figsize=(10, 5))
    plt.plot(steps, smoothed_loss_values, label='Smoothed Training Loss')
    plt.xlabel('Steps')
    plt.ylabel('Loss')
    plt.title(f'Loss Curve for Fold {fold}')
    plt.legend()
    plt.grid(True)
    plt.savefig(f"loss_curve_fold_{fold}.svg")
    plt.show()


def plot_f1_scores(f1_scores, fold):
    epochs = range(1, len(f1_scores["train"]) + 1)

    plt.figure(figsize=(10, 5))
    plt.plot(epochs, f1_scores["train"], label='Training F1 Score')
    plt.plot(epochs, f1_scores["validate"], label='Validation F1 Score')

    plt.xlabel('Epochs')
    plt.ylabel('F1 Score')
    plt.title(f'F1 Score Curve for Fold {fold}')
    plt.legend()
    plt.grid(True)
    plt.savefig(f"f1_score_curve_fold_{fold}.svg")
    plt.show()


def getFullDataLoader(X, y, batch_size):
    full_DataSet = ms.MyDataSet(input=X.reset_index(drop=True), label=y.reset_index(drop=True))
    full_DataLoader = DataLoader(dataset=full_DataSet, batch_size=batch_size, shuffle=False)
    return full_DataLoader


if __name__ == '__main__':
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    # 常用参数
    Batch_Size = 32
    test_Batch_Size = 32
    LR = 5e-4
    Epoch = 10
    K_Fold = 5
    print("Batch_Size", Batch_Size)
    print("LR", LR)
    print("Epoch", Epoch)
    print("K_Fold", K_Fold)

    file_list = ld.create_list(ld.dataDir)
    file_list.sort()
    file_list = ['cgas']

    num_gpus = torch.cuda.device_count()
    print(f"Number of GPUs available: {num_gpus}")
    device_ids = list(range(num_gpus))

    writer = SummaryWriter(log_dir='runs/experiment_1')

    for path in file_list:
        all_data = pd.read_csv(ld.dataDir + path + '/all_data.txt', sep='\t')
        print("数据列名：", all_data.columns)
        X = all_data.iloc[:, 0]  # 序列
        y = all_data.iloc[:, 1]  # 标签

        kf = StratifiedKFold(n_splits=K_Fold, shuffle=True, random_state=1)
        fold = 1
        roc_total = 0
        pr_total = 0
        F1_total = 0
        fold_results = []   # 用于绘制 ROC & PR
        metrics_list = []   # 用于保存每个 fold 的评估指标

        for train_index, validate_index in kf.split(X, y):
            train_DataLoader, validate_DataLoader, validate_Full_DataLoader = getDataSet(train_index, validate_index)

            model = BCL_Network().to(device)
            if num_gpus > 1:
                model = nn.DataParallel(model, device_ids=device_ids)

            optimizer = torch.optim.Adam(model.parameters(), lr=LR)
            scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', factor=0.5, patience=3)
            loss_func = nn.BCELoss()

            best_model_name, f1_scores, loss_values = train(train_DataLoader, validate_DataLoader, path, fold, writer)
            plot_loss_curve(loss_values, fold)
            plot_f1_scores(f1_scores, fold)

            # 🌟 收集评估指标
            ROC, PR, f1, acc, prec, rec, y_true, y_prob = evaluate(validate_Full_DataLoader, path, fold, best_model_name)

            roc_total += ROC
            pr_total += PR
            F1_total += f1

            # 🌟 收集用于绘制曲线的预测信息
            fold_results.append({
                "fold": fold,
                "y_true": y_true,
                "y_prob": y_prob
            })

            # 🌟 收集当前 fold 的指标用于保存
            metrics_list.append({
                "Fold": fold,
                "F1": f1,
                "Accuracy": acc,
                "Precision": prec,
                "Recall": rec
            })

            # 可选：保存每个 fold 的验证集预测结果
            validate(validate_DataLoader, Epoch, save_csv=True, csv_path=f"validation_predictions_fold_{fold}.csv")

            fold += 1

        # 🌟 计算平均指标
        roc_average = roc_total / K_Fold
        pr_average = pr_total / K_Fold
        f1_average = F1_total / K_Fold

        # 🌟 绘制整体 ROC 和 PR 曲线（SVG）
        plot_roc_all_folds(fold_results, path)
        # 🌟 保存 metrics 到 CSV（附加平均值）
        metrics_df = pd.DataFrame(metrics_list)
        avg_row = {
            "Fold": "Average",
            "F1": metrics_df["F1"].mean(),
            "Accuracy": metrics_df["Accuracy"].mean(),
            "Precision": metrics_df["Precision"].mean(),
            "Recall": metrics_df["Recall"].mean()
        }
        metrics_df = metrics_df.append(avg_row, ignore_index=True)
        metrics_df = metrics_df.round(4)

        metrics_out_path = os.path.join("figure", path, "metrics_per_fold.csv")
        os.makedirs(os.path.dirname(metrics_out_path), exist_ok=True)
        metrics_df.to_csv(metrics_out_path, index=False)
        print(f"✅ Saved metrics CSV to: {metrics_out_path}")

        # 🌟 控制台输出
        print(path)
        print("平均 ROC:{:.4f} | PR:{:.4f} | F1:{:.4f}".format(roc_average, pr_average, f1_average))
        print("#################################")

    # ✅ 使用最后一个最佳模型在测试集上进行预测
    for path in file_list:
        all_data = pd.read_csv(ld.dataDir + path + '/test.txt', sep='\t')
        print("数据列名：", all_data.columns)
        X = all_data.iloc[:, 0]
        y = all_data.iloc[:, 1]
        full_DataLoader = getFullDataLoader(X, y, test_Batch_Size)
        predictions, labels = predict(best_model_name, full_DataLoader)

        predictions_df = pd.DataFrame({"True Label": labels, "Predicted Label": predictions})
        predictions_df.to_csv(f"test_predictions_{path}.csv", index=False)
        print(f"✅ Test predictions saved to test_predictions_{path}.csv")

    writer.close()
