import os
import csv
import json
import time
import math
from itertools import combinations
from alignment_classic import *
from alignment_free import *


def read_fasta(file_path: str) -> str:
    with open(file_path, "r") as f:
        linee = [line.strip() for line in f if not line.startswith(">")]
    return "".join(linee)


def get_fasta_header(file_path: str) -> str:
    with open(file_path, mode="r", encoding="utf-8") as f:
        first_row = f.readline().strip()
        if first_row.startswith(">"):
            return first_row[1:]
    return "Header non trovato"


def clean_alignment(alignment_list: list, h_index: int) -> str:
    """Rimuove i None residui dall'allocazione iniziale e restituisce la stringa finale."""
    valid_chars = alignment_list[h_index:]
    return "".join([char for char in valid_chars if char is not None])


def run_pipeline(k: int = 3):
    def clean_csv_val(v):
        if isinstance(v, float):
            if math.isinf(v) or math.isnan(v):
                return str(v)
            return f"{v:.6f}"
        return v

    def clean_json_val(v):
        if isinstance(v, float) and (math.isinf(v) or math.isnan(v)):
            return str(v)
        return v

    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    CONFIG_PATH = os.path.join(BASE_DIR, "data", "pipeline_config.csv")
    OUTPUT_CSV = os.path.join(BASE_DIR, "data", "output", "allineament_results.csv")
    OUTPUT_JSON = os.path.join(BASE_DIR, "data", "output", "allineament_details.json")

    if not os.path.exists(CONFIG_PATH):
        print(f"Errore: Manca il file Manifest in {CONFIG_PATH}")
        print("Devi creare il file 'data/pipeline_config.csv' prima di avviare.")
        return

    csv_results = []
    json_results = []

    # Lettura delle famiglie dal manifest.csv
    with open(CONFIG_PATH, mode="r", newline="", encoding="utf-8") as infile:
        reader = csv.DictReader(infile)

        for riga in reader:
            group = riga["group"]
            target_dir = os.path.join(BASE_DIR, riga["directory"])

            if not os.path.isdir(target_dir):
                print(f"Saltato group {group}: La cartella {target_dir} non esiste.")
                continue

            file_genici = [f for f in os.listdir(target_dir) if f.endswith(('.fasta', '.fa'))]

            if len(file_genici) < 2:
                print(f"Avviso: Elementi insufficienti nella cartella {group}.")
                continue

            print(f"Rilevato group '{group.upper()}' con {len(file_genici)} sequenze. Generazione combinazioni...")

            for file_A, file_B in combinations(file_genici, 2):
                path_a = os.path.join(target_dir, file_A)
                path_b = os.path.join(target_dir, file_B)

                descrizione_A = get_fasta_header(path_a)
                descrizione_B = get_fasta_header(path_b)

                seq_x = read_fasta(path_a)
                seq_y = read_fasta(path_b)

                n = len(seq_x)
                m = len(seq_y)

                # Edit Distance
                ALX_ed = [None] * (n + m)
                ALY_ed = [None] * (n + m)

                start_ed = time.perf_counter()
                M_ed = edit_distance(seq_x, seq_y)
                h_index_ed = allinea(seq_x, seq_y, M_ed, ALX_ed, ALY_ed)
                end_ed = time.perf_counter()

                time_ed = end_ed - start_ed
                score_ed = M_ed[n][m]

                aligned_x_ed = clean_alignment(ALX_ed, h_index_ed)
                aligned_y_ed = clean_alignment(ALY_ed, h_index_ed)

                # Needleman-Wunsch
                ALX_nw = [None] * (n + m)
                ALY_nw = [None] * (n + m)

                start_nw = time.perf_counter()
                M_nw = needleman_wunsch(seq_x, seq_y, match=2, mismatch=-1, gap=-1)
                h_index_nw = allinea_needleman(seq_x, seq_y, M_nw, ALX_nw, ALY_nw, match=2, mismatch=-1, gap=-1)
                end_nw = time.perf_counter()

                time_nw = end_nw - start_nw
                score_nw = M_nw[n][m]

                aligned_x_nw = clean_alignment(ALX_nw, h_index_nw)
                aligned_y_nw = clean_alignment(ALY_nw, h_index_nw)

                # Smith-Waterman
                ALX_sw = [None] * (n + m)
                ALY_sw = [None] * (n + m)

                start_sw = time.perf_counter()
                M_sw = smith_waterman(seq_x, seq_y, match=2, mismatch=-1, gap=-1)
                h_index_sw, start_i_sw, start_j_sw = allinea_smith(seq_x, seq_y, M_sw, ALX_sw, ALY_sw, match=2,
                                                                   mismatch=-1, gap=-1)
                end_sw = time.perf_counter()

                time_sw = end_sw - start_sw
                score_sw = max(max(riga) for riga in M_sw)

                aligned_x_sw = clean_alignment(ALX_sw, h_index_sw)
                aligned_y_sw = clean_alignment(ALY_sw, h_index_sw)

                # --- Alignment-Free Distance Calculations ---
                
                # Euclidean Distance
                start_euclidean = time.perf_counter()
                score_euclidean = euclidean_distance(seq_x, seq_y, k)
                end_euclidean = time.perf_counter()
                time_euclidean = end_euclidean - start_euclidean

                # D2 Distance
                start_d2 = time.perf_counter()
                score_d2 = d2_distance(seq_x, seq_y, k)
                end_d2 = time.perf_counter()
                time_d2 = end_d2 - start_d2

                # Shannon Entropy
                start_shannon = time.perf_counter()
                score_shannon = shannon_entropy(seq_x, seq_y, k)
                end_shannon = time.perf_counter()
                time_shannon = end_shannon - start_shannon

                # Average Common Substring (ACS)
                start_acs = time.perf_counter()
                score_acs = average_common_substring(seq_x, seq_y)
                end_acs = time.perf_counter()
                time_acs = end_acs - start_acs

                # Normalized Compression Distance (NCD)
                start_ncd = time.perf_counter()
                score_ncd = normalized_compression_distance(seq_x, seq_y)
                end_ncd = time.perf_counter()
                time_ncd = end_ncd - start_ncd

                csv_results.append({
                    "group": group,
                    "seq_A": descrizione_A,
                    "seq_B": descrizione_B,
                    "len_A": n,
                    "len_B": m,
                    "edit_distance": score_ed,
                    "time_ed": f"{time_ed:.6f}",
                    "nw_score": score_nw,
                    "time_nw": f"{time_nw:.6f}",
                    "sw_score": score_sw,
                    "i_index_sw": start_i_sw,
                    "j_index_sw": start_j_sw,
                    "time_sw": f"{time_sw:.6f}",
                    "euclidean_dist": clean_csv_val(score_euclidean),
                    "time_euclidean": f"{time_euclidean:.6f}",
                    "d2_dist": clean_csv_val(score_d2),
                    "time_d2": f"{time_d2:.6f}",
                    "shannon_dist": clean_csv_val(score_shannon),
                    "time_shannon": f"{time_shannon:.6f}",
                    "acs_dist": clean_csv_val(score_acs),
                    "time_acs": f"{time_acs:.6f}",
                    "ncd_dist": clean_csv_val(score_ncd),
                    "time_ncd": f"{time_ncd:.6f}"
                })

                json_results.append({
                    "group": group,
                    "sequences": {
                        "seq_A_header": descrizione_A,
                        "seq_B_header": descrizione_B,
                        "len_A": n,
                        "len_B": m
                    },
                    "edit_distance_analysis": {
                        "score": score_ed,
                        "time_seconds": time_ed,
                        "h_index": h_index_ed,
                        "aligned_A": aligned_x_ed,
                        "aligned_B": aligned_y_ed
                    },
                    "needleman_wunsch_analysis": {
                        "score": score_nw,
                        "time_seconds": time_nw,
                        "h_index": h_index_nw,
                        "aligned_A": aligned_x_nw,
                        "aligned_B": aligned_y_nw
                    },
                    "smith_waterman_analysis": {
                        "score": score_sw,
                        "time_seconds": time_sw,
                        "h_index": h_index_sw,
                        "start_coordinate_i": start_i_sw,
                        "start_coordinate_j": start_j_sw,
                        "aligned_A": aligned_x_sw,
                        "aligned_B": aligned_y_sw
                    },
                    "alignment_free_analysis": {
                        "k_value": k,
                        "euclidean_distance": {
                            "score": clean_json_val(score_euclidean),
                            "time_seconds": time_euclidean
                        },
                        "d2_distance": {
                            "score": clean_json_val(score_d2),
                            "time_seconds": time_d2
                        },
                        "shannon_entropy": {
                            "score": clean_json_val(score_shannon),
                            "time_seconds": time_shannon
                        },
                        "average_common_substring": {
                            "score": clean_json_val(score_acs),
                            "time_seconds": time_acs
                        },
                        "normalized_compression_distance": {
                            "score": clean_json_val(score_ncd),
                            "time_seconds": time_ncd
                        }
                    }
                })

    os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)

    intestazioni_csv = [
        "group", "seq_A", "seq_B", "len_A", "len_B",
        "edit_distance", "time_ed",
        "nw_score", "time_nw",
        "sw_score", "i_index_sw", "j_index_sw", "time_sw",
        "euclidean_dist", "time_euclidean",
        "d2_dist", "time_d2",
        "shannon_dist", "time_shannon",
        "acs_dist", "time_acs",
        "ncd_dist", "time_ncd"
    ]
    with open(OUTPUT_CSV, mode="w", newline="", encoding="utf-8") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=intestazioni_csv)
        writer.writeheader()
        writer.writerows(csv_results)

    with open(OUTPUT_JSON, mode="w", encoding="utf-8") as json_file:
        json.dump(json_results, json_file, indent=4, ensure_ascii=False)

    print(f"Processo completato con successo!")
    print(f" -> Tabelle statistiche (CSV): {OUTPUT_CSV}")
    print(f" -> Dettagli allineamenti (JSON): {OUTPUT_JSON}")


if __name__ == "__main__":
    run_pipeline()
