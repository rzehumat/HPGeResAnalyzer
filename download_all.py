import urllib
import os
import pandas as pd
import re
import decimal

from math import pow
from pprint import pprint
from bs4 import BeautifulSoup
from termcolor import colored

Z_MIN = 51
Z_MAX = 112
TIME_CONVERSION = {
    "us": 0.000001, 
    "ms": 0.001, 
    "s": 1, 
    "m": 60, 
    "h": 3600, 
    "d": 86400, 
    "y": 31556952
}

def download_isotopes_list(Z):
    request_url = f"http://nucleardata.nuclear.lu.se/toi/listnuc.asp?sql=&Z={Z}"
    urllib.request.urlretrieve(request_url, f"downloads/find_isotopes/z_{Z}.html")

def download_range(z_min, z_max):
    for Z in range(z_min, z_max + 1):
        download_isotopes_list(Z)

def parse_isotopes_one(Z):
    html_path = f"downloads/find_isotopes/z_{Z}.html"
    isotopes_lst_html = open(html_path, "r")
    bs = BeautifulSoup(isotopes_lst_html.read(), 'lxml')
    table = bs.find_all("table")[0]
    nuclide_lst = table.find_all('th')[9:]
    abbr = str(nuclide_lst[0].find('a')).split('</sup>')[1][:-4]
    out_file_path = f"downloads/find_isotopes_parsed/{Z}_{abbr}.txt"
    out_file = open(out_file_path, "w")
    out_file.write(f"{Z}\n")
    out_file.write(f"{abbr}\n")
    for nuclide in nuclide_lst:
        out_file.write(f"{nuclide.find('sup').get_text()}\n")

    out_file.close()

def parse_isotopes_range(z_min, z_max):
    for Z in range(z_min, z_max + 1):
        parse_isotopes_one(Z)

def process_A(A):
    A = A.strip()
    if str(A)[-1] == "m":
        str_A = "0" + str(int(A[:-1]) + 300)
    elif str(A)[-2:] == "m2":
        str_A = "0" + str(int(A[:-2]) + 600)
    elif str(A)[-2:] == "m3":
        str_A = int(A[:-2]) + 900
    elif str(A)[-2:] == "m4":
        str_A = int(A[:-2]) + 1200
    elif int(A) < 10:
        str_A = '000' + A
    elif int(A) < 100:
        str_A = '00' + A
    else:
        str_A = '0' + A
    
    return A, str_A

def no_gammas_warning(A, Z, str_A, element):
    print(colored(f"Seems like there are no gamma-lines known for isotope {A}{element}.", 'yellow'))
    print(colored("Check yellow pages for reference.", 'yellow'))
    print(colored(f"http://nucleardata.nuclear.lu.se/toi/nuclide.asp?iZA={Z}{str_A}", 'yellow'))

def download_all_isotopes(Z):
    down_dir = "downloads/find_isotopes_parsed"
    files_lst = os.listdir(down_dir)
    file_name = [s for s in files_lst if str(Z) == s.split('_')[0]]
    isotopes_lst_file = open(f"{down_dir}/{file_name[0]}", "r")
    lines = isotopes_lst_file.readlines()
    abbr = lines[1].strip()
    A_lst = lines[2:]
    for A in A_lst:
        print(f"Z is {Z}")
        print(f"A is {A}")
        
        A, str_A = process_A(A)
        url = f"http://nucleardata.nuclear.lu.se/toi/nuclide.asp?iZA={Z}{str_A}"
        urllib.request.urlretrieve(url, f"downloads/isotopes_html/{A}_{abbr}_{Z}.html")

def download_all_elements(z_min, z_max):
    for Z in range(z_min, z_max + 1):
        download_all_isotopes(Z)

def clean_number(str_num):
    for char in ["<", ">", "(calc)", "*", "~", "(syst)"]:
        str_num = str(str_num).replace(char, "")
    return str_num

def extract_Igamma(A, element, Z):
    html_file = open(f"downloads/isotopes_html/{A}_{element}_{Z}.html", "r")
    soup = BeautifulSoup(html_file.read(), 'lxml')

    try:
        gammas_table = soup.find_all("table")[4]

        if any(x in gammas_table.find("font").get_text(strip=True) for x in ["Betas", "X-rays"]):
            A, str_A = process_A(A)
            no_gammas_warning(A, Z, str_A, element)
            
            return 1

        gammas_rows = gammas_table.find_all('tr')[3:-1]
    except:
        A, str_A = process_A(A)
        no_gammas_warning(A, Z, str_A, element)

        return 1
    
    energy = []
    sigm_energy = []
    i = []
    sigm_i = []   

    for row in gammas_rows:
        cells = row.find_all('td')
        
        e_cell = cells[0]
        i_cell = cells[1]
        
        if "<i>" in str(i_cell):
            sigm_ig_val = i_cell.find("i").get_text(strip=True)
            sigm_ig_val = clean_number(sigm_ig_val)
          
            ig_val = str(i_cell).split("<i>")[0].split("<td>")[1].strip()
            ig_val = clean_number(ig_val)

            ig_decimals = decimal.Decimal(ig_val)
            sigm_ig_val = float(sigm_ig_val) * pow(10, ig_decimals.as_tuple().exponent)
            ig_val = float(ig_val)

        else:
            sigm_ig_val = float("NaN")
            ig_val = i_cell.get_text(strip=True)
            ig_val = clean_number(ig_val)
            if ig_val == "":
                ig_val = float("NaN")
            else:
                ig_val = float(ig_val)
        
        if "<i>" in str(e_cell):
            sigm_e_val = e_cell.find("i").get_text(strip=True)
            sigm_e_val = clean_number(sigm_e_val)

            e_val = str(e_cell).split("<i>")[0].split("<td>")[1].strip()
            e_val = clean_number(e_val)

            e_decimals = decimal.Decimal(e_val)
            sigm_e_val = float(sigm_e_val) * pow(10, e_decimals.as_tuple().exponent)
            e_val = float(e_val)

        else:
            sigm_e_val = float("NaN")
            e_val = e_cell.get_text(strip=True)
            e_val = clean_number(e_val)

        energy.append(e_val)
        sigm_energy.append(sigm_e_val)
        i.append(ig_val)  
        sigm_i.append(sigm_ig_val)


    df_dict = {
        "E_tab": energy,
        "sigm_E": sigm_energy, 
        "Ig": i,
        "sigm_Ig": sigm_i
        }
    df = pd.DataFrame(df_dict)
    A = str(A).strip()
    df_name = f'downloads/ig_db/{A}{element}.csv'
    df["Ig"] = 0.01*df["Ig"]
    df["sigm_Ig"] = 0.01*df["sigm_Ig"]
    df.to_csv(df_name)
   
    print(f"Ig extracted from file 'downloads/{A}{element}.html' into '{df_name}'.")
    return 0

def extract_all_elements(z_min, z_max):
    for Z in range(z_min, z_max + 1):
        extract_element(Z)

def extract_element(Z):
    html_lst = os.listdir("downloads/isotopes_html")
    element_files = [f for f in html_lst if str(Z) == f.split('_')[-1].split('.')[0]]
    for isotope_file in element_files:
        A, element, Z = (isotope_file.split('.')[0]).split('_')
        print(f"Extracting {A}{element}")
        extract_Igamma(A, element, Z)
        extract_info(A, element, Z)

def extract_info(A, element, Z):
    html_file = open(f"downloads/isotopes_html/{A}_{element}_{Z}.html", "r")
    soup = BeautifulSoup(html_file.read(), 'lxml')
    table = soup.find_all("table")[0]

    info_rows = table.find_all("tr")[6:16]
    info_df = {}
    for row in info_rows:
        try:
            key = (row.find_all("th")[0]).get_text(strip=True)
        except:
            break
        try:
            val = row.find_all("td")[0]
        except:
            continue
        if val.find("i"):
            sigm = val.find("i").get_text(strip=True)
            val = str(val.get_text(strip=True))[:-len(sigm)]
            info_df[f"sigm_{key}"] = sigm
        else:
            val = val.get_text(strip=True)

        info_df[key] = val
    
    info_df.pop('', None)
    info_df = pd.DataFrame(info_df, index=[0])

    columns = info_df.columns.tolist()
    for i in range(len(columns)):
        columns[i] = columns[i][:-1]
        columns[i] = (columns[i]).replace(u'\xa0', u' ')

    info_df.columns = columns
    if "Literature cut-off date" in columns:
        info_df["Literature cut-off date"] = pd.to_datetime(info_df["Literature cut-off date"])

    if "Prod. mode" in columns:
        for mode in re.findall('[A-Z][^A-Z]*', str(info_df["Prod. mode"][0])):
            mode = mode.replace(u'\xa0', u' ')
            info_df[f"Prod_mode_{mode}"] = True
    
    for var in ["Sn(keV)", "Sp(keV)"]:
        if var in columns:
            orig_num = info_df[var][0]
            if info_df[var][0] == "":
                continue

            info_df[var] = info_df[var].astype(float)

            sigma = f"sigm_{var}"
            if sigma in columns and info_df[sigma][0] != "sy":
                g = decimal.Decimal(orig_num)
                info_df[sigma] = int(info_df[sigma]) * pow(10, g.as_tuple().exponent)

    if info_df["Half life"][0] == "stable" or info_df["Half life"][0] == "":
        info_df["Stable"] = info_df["Half life"][0] == "stable"
        info_df = info_df.drop(columns=["Prod. mode", "Half life"], errors="ignore")
        info_df.to_csv(f"downloads/ig_db/info_{A}{element}.csv")

    else:
        info_df["Stable"] = False
        hl_val, hl_unit = info_df["Half life"][0].split()

        if not hl_unit in list(TIME_CONVERSION.keys()):
            info_df["Half-life [s]"] = None
        else:
            if hl_val[0] == ">" or hl_val[0] == "~" or hl_val[0] == "<":
                hl_val = hl_val[1:]
            print(f"hl val is {hl_val}")
            d = decimal.Decimal(hl_val)
            hl_val = float(hl_val)

            info_df["Half-life [s]"] = hl_val * TIME_CONVERSION[hl_unit]
            if "sigm_Half life" in columns:
                if "+" in str(info_df["sigm_Half life"][0]):
                    info_df["sigm_Half life"][0] = ((str(info_df["sigm_Half life"][0]).split('+')[1]).split('-')[0]).strip()

                info_df["sigm_Half-life [s]"] = int(info_df["sigm_Half life"][0]) * pow(10, d.as_tuple().exponent) * TIME_CONVERSION[hl_unit]
        
        info_df = info_df.drop(columns=["Prod. mode", "Half life"], errors="ignore")
        if "sigm_Half life" in columns:
            info_df = info_df.drop(columns=["sigm_Half life"])
        if "sigm_Sp(keV)" in columns:
            info_df = info_df.drop(columns=["sigm_Sp(keV)"])
        info_df.to_csv(f"downloads/ig_db/info_{A}{element}.csv")

if __name__ == "__main__":
    download_range(Z_MIN, Z_MAX)
    parse_isotopes_range(Z_MIN,Z_MAX)
    download_all_elements(Z_MIN, Z_MAX)
    extract_all_elements(Z_MIN, Z_MAX)