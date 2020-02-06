keys = []

def getMetadata(arq):
    with open(arq, "r") as f:
        metadata = f.readline()
    metadata = [x.split('=') for x in metadata.split()]
    print(metadata)

if __name__ == "__main__":
    import sys
    getMetadata(sys.argv[1])
