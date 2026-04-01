import pandas as pd
import matplotlib.pyplot as plt

df_up = pd.read_csv("cpUp_dataTEST.csv")
df_low = pd.read_csv("cpLow_dataTEST.csv")

plt.scatter(df_low["x"], df_low["cp"], label="Lower Surface", color="blue")
plt.scatter(df_up["x"],  df_up["cp"],  label="Upper Surface", color="red")
plt.xlabel("Distance from Leading Edge")
plt.ylabel("Cp")
plt.title("Pressure Coefficient Distribution")
plt.gca().invert_yaxis()
plt.legend()
plt.savefig("cp_plot.png")
plt.show()