import requests
import time
import sqlite3
import pandas as pd


def lookup_data(url):
    web_page = requests.get(url)
    time.sleep(3)
    return web_page.text


def separate_data(info):
    info = info.replace("INFINITE", "inf")
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
    compound_formula = compound_id.split(" (")[1].split(")	")[0]
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

    df = pd.DataFrame(
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
    return df, compound_name, compound_formula


def store_database(c_name, c_formula, df1):
    db_conn = sqlite3.connect("NIST_JANAF_Database.sqlite")
    cur = db_conn.cursor()
    # cur.execute("CREATE TABLE sample (formula VARCHAR, T_K VARCHAR, Cp VARCHAR, S VARCHAR, G_H_Tr_T VARCHAR,"
    #            "H_H_Tr VARCHAR, Df_H VARCHAR, Df_G VARCHAR, log_K VARCHAR)")
    # conn.commit()
    df1.to_sql(cur)
    # cur.execute("INSERT INTO sample (formula, T_K, Cp, S, G_H_Tr_T, H_H_Tr, Df_H, Df_G, log_K)"
    #            "VALUES (?, ?, ?, ?, ?, ?, ?, ?), (c_formula, df1['T_K'], df1['Cp'], df1['S'],"
    #            "df1['G_H_Tr_T'], df1['H_H_Tr'], df1['Df_H'], df1['Df_G'], df1['log_K']))")
    conn.commit()


sample_data = """Niobium Carbide (NbC0.98)	C0.98Nb1(cr)
T(K)	Cp	S	-[G-H(Tr)]/T	H-H(Tr)	delta-f H	delta-f G	log Kf
0	0.	0.	INFINITE	-5.422	-138.061	-138.061	INFINITE
100	14.426	7.983	56.725	-4.874	-138.352	-137.801	71.980
200	27.476	22.238	35.932	-2.739	-138.635	-137.103	35.808
298.15	36.233	34.966	34.966	0.	-138.909	-136.785	23.964
300	36.380	35.190	34.966	0.067	-138.903	-136.772	23.814
400	41.882	46.463	36.466	3.999	-138.481	-136.120	17.776
500	45.145	56.179	39.461	8.359	-137.986	-135.587	14.165
600	47.363	64.616	42.966	12.990	-137.514	-135.153	11.766
700	48.911	72.038	46.600	17.806	-137.091	-134.793	10.058
800	50.124	78.651	50.200	22.760	-136.718	-134.491	8.781
900	51.087	84.611	53.698	27.822	-136.389	-134.233	7.791
1000	51.923	90.038	57.065	32.973	-136.097	-134.009	7.000
1100	52.635	95.021	60.292	38.202	-135.839	-133.813	6.354
1200	53.304	99.630	63.380	43.499	-135.609	-133.640	5.817
1300	53.932	103.921	66.336	48.861	-135.398	-133.484	5.363
1400	54.504	107.939	69.165	54.283	-135.207	-133.344	4.975
1500	55.061	111.718	71.878	59.761	-135.038	-133.217	4.639
1600	55.605	115.289	74.480	65.295	-134.894	-133.101	4.345
1700	56.134	118.676	76.981	70.882	-134.781	-132.992	4.086
1800	56.651	121.900	79.388	76.521	-134.705	-132.889	3.856
1900	57.153	124.976	81.707	82.211	-134.671	-132.790	3.651
2000	57.656	127.920	83.944	87.952	-134.683	-132.691	3.466
2100	58.158	130.746	86.106	93.743	-134.745	-132.590	3.298
2200	58.660	133.463	88.197	99.583	-134.863	-132.484	3.146
2300	59.162	136.081	90.223	105.475	-135.042	-132.373	3.006
2400	59.664	138.610	92.187	111.416	-135.293	-132.251	2.878
2500	60.166	141.056	94.093	117.407	-135.630	-132.118	2.760
2600	60.668	143.425	95.945	123.449	-136.076	-131.969	2.651
2700	61.170	145.724	97.746	129.541	-136.669	-131.801	2.550
2800	61.672	147.958	99.500	135.683	-163.822	-131.122	2.446
2900	62.174	150.131	101.208	141.875	-163.565	-129.959	2.341
3000	62.676	152.247	102.874	148.118	-163.272	-128.805	2.243"""


conn = sqlite3.connect("NIST_JANAF_Database.sqlite")
conn.close()

for i in range(0, 200):
    num_id = str(i).zfill(3)
    URL = f"https://janaf.nist.gov/tables/C-{num_id}.txt"
    # raw_table = lookup_data(URL)
    raw_table = sample_data
    if "404 Not Found" in raw_table:
        print("Loop done")
        break
    df, compound_name, compound_formula = separate_data(raw_table)
    # store_database(df, compound_name, compound_formula)
    print(df)
