import pandas as pd
import matplotlib.pyplot as plt

def main():
    total = pd.read_csv("//wsl.localhost/Ubuntu/home/angel/total_method.csv")
    conv  = pd.read_csv("//wsl.localhost/Ubuntu/home/angel/convection_method.csv")

    plt.plot(total["mf"], total["dGC_total"], label="Total")
    plt.plot(conv["mf"], conv["dGC_convection"], label="Convection")

    # Etiquetas
    plt.xlabel(r"$m_{\varphi}$")
    plt.ylabel(r"$\Gamma(\varphi \rightarrow K^+ K^-)$")

    plt.legend()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()