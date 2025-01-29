import os, sys

from DataLoader import MathProblemsDataset

def main():
    # Get data
    ds = MathProblemsDataset(dataset_name="AI-MO/NuminaMath-TIR", partition_name='train')

    dataset = ds.get_data()
    problems = dataset['problem']

    print(problems[:4])


if __name__ == "__main__":
    print('Running ...')
    main()