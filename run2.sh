#!/bin/bash

# --- CONFIGURACIÓN ---
# Tamaños de grid a probar para el análisis de escalabilidad
GRID_SIZES_A_PROBAR="2048 4096 8192 16384" 
PROCESOS_A_PROBAR="1 2 4 8"
THREADS_A_PROBAR="1 2 4 8"
OUTPUT_CSV="resultados_escalabilidad_completo.csv"

# --- COMPILACIÓN (Solo una vez) ---
echo "--- Iniciando compilación de todos los códigos ---"
g++ Secuencial/marching_square.cpp -o Secuencial/a.out || { echo "Error compilando Secuencial"; exit 1; }
g++-14 -fopenmp Compartido/marching_squares.cpp -o Compartido/a.out || { echo "Error compilando Compartido"; exit 1; }
mpic++ Distribuido/marching_squares.cpp -o Distribuido/a.out || { echo "Error compilando Distribuido"; exit 1; }
mpic++ -fopenmp Hibrido/marching_squares.cpp -o Hibrido/a.out || { echo "Error compilando Híbrido"; exit 1; }
echo "--- Compilación finalizada exitosamente ---"
echo ""

# --- EJECUCIÓN Y RECOLECCIÓN DE DATOS ---
echo "--- Iniciando pruebas de escalabilidad ---"
echo "Resultados se guardarán en: $OUTPUT_CSV"

# Crear el archivo CSV y escribir los nuevos encabezados
echo "grid_size,paradigma,procesos,threads_por_proceso,total_threads,tiempo_ms" > $OUTPUT_CSV

# Bucle externo para cada tamaño de grid
for gs in $GRID_SIZES_A_PROBAR; do
    echo ""
    echo "*****************************************************"
    echo "--- PROBANDO GRID SIZE: ${gs}x${gs} ---"
    echo "*****************************************************"

    # 1. Prueba Secuencial (nuestro baseline para cada grid_size)
    echo "Ejecutando: Secuencial (Grid: ${gs}x${gs})"
    TIME_MS_SEQ=$(./Secuencial/a.out $gs $gs | grep "FINAL_TIME_MS:" | cut -d':' -f2)
    echo "$gs,Secuencial,1,1,1,$TIME_MS_SEQ" >> $OUTPUT_CSV

    # 2. Pruebas Compartido (OpenMP)
    for h in $THREADS_A_PROBAR; do
        if [ $h -eq 1 ]; then
            echo "$gs,Compartido,1,1,1,$TIME_MS_SEQ" >> $OUTPUT_CSV
            continue
        fi
        echo "Ejecutando: Compartido (Grid: ${gs}x${gs}, Threads: $h)"
        TIME_MS=$(./Compartido/a.out $gs $gs $h | grep "FINAL_TIME_MS:" | cut -d':' -f2)
        echo "$gs,Compartido,1,$h,$h,$TIME_MS" >> $OUTPUT_CSV
    done

    # 3. Pruebas Distribuido (MPI)
    for p in $PROCESOS_A_PROBAR; do
        if [ $p -eq 1 ]; then
            echo "$gs,Distribuido,1,1,1,$TIME_MS_SEQ" >> $OUTPUT_CSV
            continue
        fi
        echo "Ejecutando: Distribuido (Grid: ${gs}x${gs}, Procesos: $p)"
        TIME_MS=$(mpirun -np $p ./Distribuido/a.out $gs $gs | grep "Tiempo de Cómputo:" | awk '{print $4}')
        echo "$gs,Distribuido,$p,1,$p,$TIME_MS" >> $OUTPUT_CSV
    done

    # 4. Pruebas Híbrido (MPI + OpenMP)
    for p in $PROCESOS_A_PROBAR; do
        for h in $THREADS_A_PROBAR; do
            TOTAL_THREADS=$((p * h))
            if [ $p -eq 1 ] && [ $h -eq 1 ]; then
                echo "$gs,Hibrido,1,1,1,$TIME_MS_SEQ" >> $OUTPUT_CSV
                continue
            fi
            echo "Ejecutando: Hibrido (Grid: ${gs}x${gs}, Procesos: $p, Threads: $h)"
            TIME_MS=$(mpirun -np $p ./Hibrido/a.out $gs $gs $h | grep "Tiempo de Cómputo:" | awk '{print $4}')
            echo "$gs,Hibrido,$p,$h,$TOTAL_THREADS,$TIME_MS" >> $OUTPUT_CSV
        done
    done
done

echo ""
echo "--- Script finalizado. ¡Pruebas de escalabilidad completas! ---"
