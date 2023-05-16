from quinoa import *

def main():
    rho = state.DensityMatrix.vacuum()
    for c, ketbra in rho.ketbras():
        print(c, ketbra)

if __name__ == "__main__":
    main()
