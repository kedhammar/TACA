import pandas as pd


def oops():
    pd.DataFrame({"a": [1, 2, 3],    "b": [4, 5, 6]}).to_csv("test.csv")
