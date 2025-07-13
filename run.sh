#!/bin/bash

# --- CONFIGURACIÓN ---
GRID_SIZE=4096      # Tamaño de la grilla para todas las pruebas (ej. 1024x1024)
MAX_PROCESOS=4     # Número máximo de procesos MPI a probar
MAX_THREADS=4       # Número máximo de hilos OpenMP a probar
OUTPUT_CSV="resultados_sin_opt.csv"

# --- COMPILACIÓN (SIN -O3) ---
echo "--- Iniciando compilación de todos los códigos (SIN optimización -O3) ---"

# Compilar Secuencial
g++ Secuencial/marching_square.cpp -o Secuencial/a.out || { echo "Error compilando Secuencial"; exit 1; }

# Compilar Compartido (OpenMP)
g++ -fopenmp Compartido/marching_squares.cpp -o Compartido/a.out || { echo "Error compilando Compartido"; exit 1; }

# Compilar Distribuido (MPI)
mpic++ Distribuido/marching_squares.cpp -o Distribuido/a.out || { echo "Error compilando Distribuido"; exit 1; }

# Compilar Híbrido (MPI + OpenMP)
mpic++ -fopenmp Hibrido/marching_squares.cpp -o Hibrido/a.out || { echo "Error compilando Híbrido"; exit 1; }

echo "--- Compilación finalizada exitosamente ---"
echo ""


# --- EJECUCIÓN Y RECOLECCIÓN DE DATOS ---
echo "--- Iniciando pruebas de rendimiento ---"
echo "Resultados se guardarán en: $OUTPUT_CSV"

# Crear el archivo CSV y escribir los encabezados
echo "paradigma,procesos,threads_por_proceso,total_threads,tiempo_ms" > $OUTPUT_CSV

# 1. Prueba Secuencial (nuestro baseline para 1,1)
echo "Ejecutando: Secuencial (1 Proceso, 1 Thread)"
TIME_MS_SEQ=$(./Secuencial/a.out $GRID_SIZE $GRID_SIZE | grep "FINAL_TIME_MS:" | cut -d':' -f2)
echo "Secuencial,1,1,1,$TIME_MS_SEQ" >> $OUTPUT_CSV


# 2. Pruebas Compartido (OpenMP)
for h in $(seq 1 $MAX_THREADS); do
    if [ $h -eq 1 ]; then
        echo "Compartido,1,1,1,$TIME_MS_SEQ" >> $OUTPUT_CSV
        continue
    fi
    echo "Ejecutando: Compartido (1 Proceso, $h Threads)"
    TIME_MS=$(./Compartido/a.out $GRID_SIZE $GRID_SIZE $h | grep "FINAL_TIME_MS:" | cut -d':' -f2)
    echo "Compartido,1,$h,$h,$TIME_MS" >> $OUTPUT_CSV
done


# 3. Pruebas Distribuido (MPI)
for p in $(seq 1 $MAX_PROCESOS); do
     if [ $p -eq 1 ]; then
        echo "Distribuido,1,1,1,$TIME_MS_SEQ" >> $OUTPUT_CSV
        continue
    fi
    echo "Ejecutando: Distribuido ($p Procesos, 1 Thread)"
    TIME_MS=$(mpirun -np $p ./Distribuido/a.out $GRID_SIZE $GRID_SIZE | grep "Tiempo de Cómputo:" | awk '{print $4}')
    echo "Distribuido,$p,1,$p,$TIME_MS" >> $OUTPUT_CSV
done


# 4. Pruebas Híbrido (MPI + OpenMP)
for p in $(seq 1 $MAX_PROCESOS); do
    for h in $(seq 1 $MAX_THREADS); do
        TOTAL_THREADS=$((p * h))
        if [ $p -eq 1 ] && [ $h -eq 1 ]; then
            echo "Hibrido,1,1,1,$TIME_MS_SEQ" >> $OUTPUT_CSV
            continue
        fi

        echo "Ejecutando: Hibrido ($p Procesos, $h Threads)"
        TIME_MS=$(mpirun -np $p ./Hibrido/a.out $GRID_SIZE $GRID_SIZE $h | grep "Tiempo de Cómputo:" | awk '{print $4}')
        echo "Hibrido,$p,$h,$TOTAL_THREADS,$TIME_MS" >> $OUTPUT_CSV
    done
done

echo ""
echo "--- Script finalizado. ¡Pruebas completas! ---"
