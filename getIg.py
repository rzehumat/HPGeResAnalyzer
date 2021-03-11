import pandas as pd
import json
import numpy as np
import os
import glob
import re
import urllib
import shutil

from bs4 import BeautifulSoup
from pathlib import Path
from sys import argv
from termcolor import colored

def append_Igamma_dir(parsed_dir):
    if not os.path.isdir(parsed_dir):
        raise Exception(f"Directory '{parsed_dir}' not found. Consider creating it and moving parsed files there.")

    for parsed_file in glob.iglob(f"{parsed_dir}/*.csv"):
        append_Igamma(parsed_file)

def append_Igamma(file):
    A, element = re.split(r'(\d+)(\w+)', file)[-3:-1]
    print(f"File is {file}")
    print(A)
    element = element.split('_g')[0]
    print(element)

    Path("./ig_db").mkdir(parents=True, exist_ok=True)
    if os.path.isfile(f"ig_db/{A}{element}.csv"):
        ig_df = pd.read_csv(f"ig_db/{A}{element}.csv", header=0, index_col=0)
        print("Reading the csv with gammas.")
    else:
        Path("./downloads").mkdir(parents=True, exist_ok=True)
        if os.path.isfile(f"downloads/{A}{element}.html"):
            err = extract_Igamma(A, element)
            if err > 0:
                shutil.copy(file, "with_Ig/")
                return "Are you sure such element has gamma lines?"
            ig_df = pd.read_csv(f"ig_db/{A}{element}.csv", header=0, index_col=0)
            print("Reading the csv with gammas.")
        else:
            download(A, element)
            err = extract_Igamma(A, element)
            if err > 0:
                shutil.copy(file, "with_Ig/")
                return "Are you sure such element has gamma lines?"
            ig_df = pd.read_csv(f"ig_db/{A}{element}.csv", header=0, index_col=0)
            print("Reading the csv with gammas.")
    
    parsed_df = pd.read_csv(file, index_col=0)
    joined_df = add_Ig(parsed_df, ig_df)
    
    Path("./with_Ig").mkdir(parents=True, exist_ok=True)
    joined_df.to_csv(f"with_Ig/{file.split('/')[-1]}")

def add_Ig(df, ig):
    a = np.empty((df.shape[0], 2))
    a[:] = np.nan
    df[["E_tab", "Ig"]] = a
    # UGLY, never for-loop in pandas
    for row_no in range(df.shape[0]):
        for ig_row_no in range(ig.shape[0]):
            if (df["Energy"] - df["FWHM"]).iloc[row_no] > ig["E_tab"].iloc[ig_row_no]:
                continue
            elif (df["Energy"] + df["FWHM"]).iloc[row_no] < ig["E_tab"].iloc[ig_row_no]:
                break
            else:
                df.loc[row_no, ["E_tab", "Ig"]] = ig.loc[ig_row_no, ["E_tab", "Ig"]]
    print("Ig added")
    return df

def getZ(element):
    element_Z ={
        "H": 1,
        "He": 2,
        "Xe": 54,
        "Eu": 63,
        "U": 92,
        "Pu": 94
    }
    return element_Z[element]

def download(A, element):
    A = int(A)
    if A < 10:
        str_A = '00' + str(A)
    elif A < 100:
        str_A = '0' + str(A)
    else:
        str_A = str(A)

    url = f"http://nucleardata.nuclear.lu.se/toi/nuclide.asp?iZA={getZ(element)}0{str_A}"
    urllib.request.urlretrieve (url, f"downloads/{A}{element}.html")
    print(f"File {A}{element}.html downloaded.")

def extract_Igamma(A, element):
    html_file = open(f"downloads/{A}{element}.html", "r")
    soup = BeautifulSoup(html_file.read(), 'lxml')
    
    try:
        gammas_table = soup.find_all("table")[4]
        gammas_rows = gammas_table.find_all('tr')[3:-1]
    except:
        A = int(A)
        if A < 10:
            str_A = '00' + str(A)
        elif A < 100:
            str_A = '0' + str(A)
        print(colored(f"Seems like there are no gamma-lines known for isotope {A}{element}.", 'red'))
        print(colored("Check yellow pages for reference.", 'yellow'))
        print(colored(f"http://nucleardata.nuclear.lu.se/toi/nuclide.asp?iZA={getZ(element)}0{str_A}", 'yellow'))
        return 1
    energy = []
    sigm_energy = []
    i = []
    sigm_i = []

    for row in gammas_rows:
        cells = row.find_all('td')
        
        e_val = cells[0].get_text(strip=True)
        i_val = cells[1].get_text(strip=True)
        try:
            ig_val = float(i_val[:-1])
            sigm_ig_val = float(i_val[-1])
        except:
            ig_val = float('NaN')
            sigm_ig_val = float('NaN')
       
        energy.append(float(e_val[:-1]))
        sigm_energy.append(int(e_val[-1]))
        i.append(ig_val)
        
        sigm_i.append(sigm_ig_val)


    df_dict = {
        "E_tab": energy,
        "sigm_E": sigm_energy, 
        "Ig": i,
        "sigm_Ig": sigm_i
        }
    df = pd.DataFrame(df_dict)
    df_name = f'ig_db/{A}{element}.csv'
    df.to_csv(df_name)
   
    print(f"Ig extracted from file 'downloads/{A}{element}.html' into '{df_name}'.")
    return 0

if __name__ == "__main__":
    folder = argv[1]
    append_Igamma_dir(folder)