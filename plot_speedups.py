import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import sys
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

INPUT_CSV_FILE = 'resultados_escalabilidad_potencias.csv'
OUTPUT_FOLDER = 'speedups'


def load_and_prepare_data(file_path: str) -> pd.DataFrame | None:
    """Carga los datos, los valida y calcula el speedup para todas las entradas."""
    try:
        df = pd.read_csv(file_path)
        logging.info(f"Datos cargados desde '{file_path}'.")
    except FileNotFoundError:
        logging.error(f"Error: No se encontró el archivo '{file_path}'.")
        return None
    except Exception as e:
        logging.error(f"Error inesperado al leer el CSV: {e}")
        return None
    
    try:
        sequential_time = df.query("paradigma == 'Secuencial'")['tiempo_ms'].iloc[0]
        df['speedup'] = sequential_time / df['tiempo_ms']
        logging.info("Cálculo de Speedup completado para todos los datos.")
        return df
    except IndexError:
        logging.error("No se encontró el tiempo de ejecución 'Secuencial' en los datos.")
        return None

def create_hybrid_plot(df_hibrido: pd.DataFrame, output_folder: str):
    """Genera y guarda la gráfica detallada para el paradigma Híbrido."""
    if df_hibrido.empty:
        logging.warning("No hay datos para el paradigma Híbrido. Se omitirá esta gráfica.")
        return

    fig, ax = plt.subplots(figsize=(12, 8))
    unique_processes = sorted(df_hibrido['procesos'].unique())
    colors = plt.cm.viridis(np.linspace(0, 1, len(unique_processes)))

    for p, color in zip(unique_processes, colors):
        subset = df_hibrido[df_hibrido['procesos'] == p].sort_values('threads_por_proceso')
        label = f'{p} {"Proceso" if p == 1 else "Procesos"}'
        ax.plot(subset['total_threads'], subset['speedup'], marker='o', linestyle='-', color=color, label=label)

    all_threads = sorted(df_hibrido['total_threads'].unique())
    ax.plot(all_threads, all_threads, linestyle='--', color='k', alpha=0.7, label='Speedup Ideal')

    ax.set_title('Análisis de Speedup Híbrido por Nº de Procesos 🚀', fontsize=18, weight='bold')
    ax.set_xlabel('Número Total de Threads', fontsize=12)
    ax.set_ylabel('Speedup (Aceleración)', fontsize=12)
    ax.set_xscale('log', base=2)
    ax.set_xticks(all_threads)
    ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
    ax.legend(title='Configuración Híbrida', fontsize=11)
    plt.tight_layout(pad=1.5)

    output_path = os.path.join(output_folder, "speedup_hibrido_detallado.png")
    plt.savefig(output_path, dpi=150)
    logging.info(f"✅ Gráfica HÍBRIDA guardada en: {output_path}")
    plt.close(fig)

def create_simple_paradigm_plot(df: pd.DataFrame, paradigma: str, output_folder: str):
    """Genera una gráfica de speedup simple para un paradigma (Compartido o Distribuido)."""
    subset = df[df['paradigma'] == paradigma].sort_values('total_threads')
    if subset.empty:
        logging.warning(f"No hay datos para el paradigma '{paradigma}'. Se omitirá esta gráfica.")
        return
        
    fig, ax = plt.subplots(figsize=(10, 7))
    color_map = {'Compartido': 'dodgerblue', 'Distribuido': 'forestgreen'}
    
    ax.plot(subset['total_threads'], subset['speedup'], marker='o', linestyle='-', label=f'Speedup {paradigma}', color=color_map.get(paradigma, 'black'))
    
    threads_in_subset = sorted(subset['total_threads'].unique())
    ax.plot(threads_in_subset, threads_in_subset, linestyle='--', color='k', alpha=0.8, label='Speedup Ideal')

    ax.set_title(f'Análisis de Speedup ({paradigma})', fontsize=16, weight='bold')
    ax.set_xlabel('Número Total de Threads', fontsize=12)
    ax.set_ylabel('Speedup (Aceleración)', fontsize=12)
    ax.set_xscale('log', base=2)
    ax.set_xticks(threads_in_subset)
    ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
    ax.legend(title='Rendimiento', fontsize=11)
    plt.tight_layout()

    output_path = os.path.join(output_folder, f"speedup_{paradigma.lower()}.png")
    plt.savefig(output_path, dpi=120)
    logging.info(f"✅ Gráfica {paradigma.upper()} guardada en: {output_path}")
    plt.close(fig)

def main():
    """Función principal para orquestar la creación de todas las gráficas."""
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    
    df = load_and_prepare_data(INPUT_CSV_FILE)
    if df is None:
        sys.exit(1)

    df_hibrido = df[df['paradigma'] == 'Hibrido']
    create_hybrid_plot(df_hibrido, OUTPUT_FOLDER)

    for paradigma in ['Compartido', 'Distribuido']:
        create_simple_paradigm_plot(df, paradigma, OUTPUT_FOLDER)

if __name__ == "__main__":
    main()