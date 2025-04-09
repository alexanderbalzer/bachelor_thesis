import argparse




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="tba")
    parser.add_argument("--input", type=str, required= True, help="Input file path")
    args = parser.parse_args()