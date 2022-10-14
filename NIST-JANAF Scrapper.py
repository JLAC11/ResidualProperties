import requests
import time
import sqlite3
import pandas as pd
from sqlalchemy import create_engine


def lookup_data(url):
    web_page = requests.get(url)
    time.sleep(3)
    return web_page.text


def separate_data(info):
    info = info.replace("INFINITE", "0")
    T_K = []
    Cp = []
    S = []
    G_H_Tr_T = []
    H_H_Tr = []
    Df_H = []
    Df_G = []
    log_Kf = []
    compound_id = str.splitlines(info)[0]
    compound_name = compound_id.split(" (")[0]
    compound_name = compound_name.replace(" ", "_")
    compound_formula = compound_id.split(" (")[1].split(")	")[0]
    # compound_name = compound_name + "___" + compound_formula
    for j in str.splitlines(info)[2:]:
        therm_state = j.split()
        T_K.append(therm_state[0])
        Cp.append(therm_state[1])
        S.append(therm_state[2])
        G_H_Tr_T.append(therm_state[3])
        H_H_Tr.append(therm_state[4])
        Df_H.append(therm_state[5])
        Df_G.append(therm_state[6])
        log_Kf.append(therm_state[7])

    df1 = pd.DataFrame(
        {
            "T_K": T_K,
            "Cp": Cp,
            "S": S,
            "G_H_Tr_T": G_H_Tr_T,
            "H_H_Tr": H_H_Tr,
            "Df_H": Df_H,
            "Df_G": Df_G,
            "log_Kf": log_Kf,
        }
    )
    return df1, compound_name


def store_database(c_name, df1):
    engine = create_engine(
        "sqlite:///T://Adal//IBERO//OTHER//Control//Scripts//NIST_JANAF_Database.sqlite",
        echo=False,
    )
    df1.to_sql(c_name, con=engine, if_exists="replace", index=False)
    engine.execute(f"SELECT * FROM {c_name}").fetchall()


for i in range(40, 500):
    try:
        num_id = str(i).zfill(3)
        URL = f"https://janaf.nist.gov/tables/C-{num_id}.txt"
        raw_table = lookup_data(URL)
        if "404 Not Found" in raw_table:
            print("Loop done")
            break
        if "TRANSITION" in raw_table:
            print(i, "NA")
            continue
        if "LIQUID" in raw_table:
            print(i, "NA")
            continue
        df, compound_name = separate_data(raw_table)
        print(i, compound_name)
        store_database(compound_name, df)
    except:
        print(i, "NA, T/E")
        continue
