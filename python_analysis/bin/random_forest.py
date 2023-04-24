import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import date
from pathlib import Path
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier


def input(infile):
    indata = pd.DataFrame(pd.read_csv(infile, sep=',', header=0))
    df = indata.dropna()
    X = df.drop(["case", "site", "Pop", "sex"], axis=1)
    y = df["sex"]
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=44)
    rf_model = RandomForestClassifier(n_estimators=50, max_features="auto", random_state=44, bootstrap=True)
    rf_model.fit(X_train, y_train)
    predictions = rf_model.predict(X_test)
    print(X_test.shape)
    print(len(predictions))

    print(y)
   

def main():
    infile = Path(sys.argv[1])
    input(infile)
    #out_prefix = Path(sys.argv[2])


if __name__ == "__main__":
    main()      
