import rptParser
import getIg
import getEpsilon

from sys import argv



rptParser.parse_RPT(argv[1])
getIg.append_Igamma_dir("parsed_reports/")
getEpsilon.add_epsilon_dir("with_Ig/", argv[2])


# import os
# import re
# import pandas as pd
# 
# from sys import argv
# from shutil import move
# from pathlib import Path
# 
# folder = argv[1]
# 
# 
# CALIBRATION_PATH = "aux_data/epsilons.csv"
# eps_df = pd.read_csv(CALIBRATION_PATH, index_col=0)
# 
# for file_name in os.listdir(folder):
#     # try:
#     #     A, element = re.split(r'(\d+)(\w+)', file_name)[-3:-1]
#     #     print(f"Isotope {A}{element} detected from filename '{file_name}'. Is that true?")
#     #     key_input = input("[Y/N]: ")
#     # except:
#     #     key_input = "N"
# 
#     # if key_input == "N":
#     #     A, element = input("Enter A and Element separated by space, e.g. '235 U', '12 C',...: ").split(' ')
#     print(f"File '{file_name}'. What is it?")
#     A, element = input("Enter A and Element separated by space, e.g. '235 U', '12 C',...: ").split(' ')
# 
#     geom = input("What geometry was it measured in? [3/30/80/120/250]: ")
#     res_df = rptParser.parse_one_RPT(os.path.join(folder, file_name))
#     print(res_df)
# 
#     Path("./ig_db").mkdir(parents=True, exist_ok=True)
#     if os.path.isfile(f"ig_db/{A}{element}.csv"):
#         ig_df = pd.read_csv(f"ig_db/{A}{element}.csv", header=0, index_col=0)
#         print("Reading the csv with gammas.")
#     else:
#         Path("./downloads").mkdir(parents=True, exist_ok=True)
#         if os.path.isfile(f"downloads/{A}{element}.html"):
#             getIg.extract_Igamma(A, element)
#             ig_df = pd.read_csv(f"ig_db/{A}{element}.csv", header=0, index_col=0)
#             print("Reading the csv with gammas.")
#         else:
#             getIg.download(A, element)
#             getIg.extract_Igamma(A, element)
#             ig_df = pd.read_csv(f"ig_db/{A}{element}.csv", header=0, index_col=0)
#             print("Reading the csv with gammas.")
#     
#     df = getIg.add_Ig(res_df, ig_df)
#     df["eps"] = df["E_tab"].apply(lambda x: getEpsilon.add_epsilon_val(x, geom, eps_df))
# 
#     # change the order of columns 
#     old_column_order = df.columns.tolist()
#     df["Geometry"] = [geom for i in range(df.shape[0])] 
#     first_cols = ["Pk", "Energy", "FWHM", "E_tab", "Area", "%err", "Ig", "eps", "Geometry"]
#     new_cols = first_cols + list(set(old_column_order) - set(first_cols))
#     df = df[new_cols]
#     # that is wrong -- it's always 30 in the RPT files...
#     df = df.drop(columns = ["Sample Geometry"])
# 
#     Path("./out").mkdir(parents=True, exist_ok=True)
#     df.to_csv(f"out/{file_name}")