from pathlib import Path
import matplotlib.pyplot as plt

SCRIPT_DIR = str(Path(__file__).parent.absolute())
CONFIG_DIR = str(Path(__file__).parent.absolute())+"/config/"

COLORS = plt.get_cmap("tab10")