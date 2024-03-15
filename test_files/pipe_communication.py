import sys


def main():
    # Read message from stdin
    message = sys.stdin .readline().strip()
    print("Received message:", message)


if __name__ == "__main__":
    main()
