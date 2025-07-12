#!/bin/bash

# ==============================================================================
# Script para automatizar pruebas con un tiempo secuencial CANÓNICO como base.
# Todos los casos de "1 núcleo" usarán el mismo tiempo medido del ejecutable secuencial.
# ==============================================================================

# --- Configuración ---
M=1000
N=1000
GRID_SIZE="${M}x${N}"
OUTPUT_FILE="data_canonical_results.txt"
CORE_SEQUENCE="1 2 3 4"

# --- Fase 1: Compilación de todos los programas ---
echo "--- ⚙️  Iniciando Fase de Compilación General ---"

g++ -std=c++17 Secuencial/marching_square.cpp -o sequential_exec
if [ $? -ne 0 ]; then echo "❌ Error al compilar Secuencial."; exit 1; fi

g++ -std=c++17 -fopenmp Compartido/marching_squares.cpp -o shared_exec
if [ $? -ne 0 ]; then echo "❌ Error al compilar Compartido."; exit 1; fi

mpic++ -std=c++17 Distribuido/marching_squares.cpp -o distributed_exec
if [ $? -ne 0 ]; then echo "❌ Error al compilar Distribuido."; exit 1; fi

mpic++ -std=c++17 -fopenmp Hibrido/marching_squares.cpp -o hybrid_exec
if [ $? -ne 0 ]; then echo "❌ Error al compilar Híbrido."; exit 1; fi

echo "--- ✅ Compilación finalizada con éxito ---"
echo ""

# --- Fase 2: Medición del Tiempo Base y Ejecución de Pruebas ---
# Crear/Limpiar archivo de resultados y añadir encabezado
echo "paradigm,processes,threads,total_cores,grid_size,time_ms" > "$OUTPUT_FILE"
echo "--- 🚀 Iniciando Pruebas de Rendimiento Canónicas ---"
echo "Resultados se guardarán en: $OUTPUT_FILE"
echo ""

# == Medir el Tiempo Secuencial Canónico UNA SOLA VEZ ==
echo "--- 🏃‍♂️ Obteniendo Tiempo Base Secuencial Canónico ---"
TIME_SEQUENTIAL_CANONICAL=$(./sequential_exec $M $N | grep "FINAL_TIME_MS:" | awk -F: '{print $2}')
echo "Tiempo base secuencial medido: $TIME_SEQUENTIAL_CANONICAL ms"
echo ""

# Escribir la primera línea en el archivo de datos para el paradigma secuencial
echo "sequential,1,1,1,$GRID_SIZE,$TIME_SEQUENTIAL_CANONICAL" >> "$OUTPUT_FILE"


# == Pruebas de Memoria Compartida (OpenMP) ==
echo "--- 👨‍👩‍👧‍👦 Ejecutando Pruebas de MEMORIA COMPARTIDA ---"
for T in $CORE_SEQUENCE; do
    if [ "$T" -eq 1 ]; then
        echo "Usando tiempo base para 1 hilo..."
        echo "shared,1,1,1,$GRID_SIZE,$TIME_SEQUENTIAL_CANONICAL" >> "$OUTPUT_FILE"
    else
        echo "Ejecutando OMP con $T hilos..."
        TIME=$(./shared_exec $M $N $T | grep "FINAL_TIME_MS:" | awk -F: '{print $2}')
        echo "shared,1,$T,$T,$GRID_SIZE,$TIME" >> "$OUTPUT_FILE"
    fi
done
echo ""

# == Pruebas de Memoria Distribuida (MPI) ==
echo "--- 🌐 Ejecutando Pruebas de MEMORIA DISTRIBUIDA ---"
for P in $CORE_SEQUENCE; do
    if [ "$P" -eq 1 ]; then
        echo "Usando tiempo base para 1 proceso..."
        echo "distributed,1,1,1,$GRID_SIZE,$TIME_SEQUENTIAL_CANONICAL" >> "$OUTPUT_FILE"
    else
        echo "Ejecutando MPI con $P procesos..."
        TIME=$(mpirun --oversubscribe -np $P ./distributed_exec $M $N | grep "Total time:" | awk '{print $3}')
        echo "distributed,$P,1,$P,$GRID_SIZE,$TIME" >> "$OUTPUT_FILE"
    fi
done
echo ""

# == Pruebas Híbridas (MPI + OpenMP) ==
echo "--- ⚛️  Ejecutando Pruebas HÍBRIDAS ---"
for P in $CORE_SEQUENCE; do
    for T in $CORE_SEQUENCE; do
        TOTAL_CORES=$((P * T))
        if [ "$P" -eq 1 ] && [ "$T" -eq 1 ]; then
            echo "Usando tiempo base para P=1, T=1..."
            echo "hybrid,1,1,1,$GRID_SIZE,$TIME_SEQUENTIAL_CANONICAL" >> "$OUTPUT_FILE"
        else
            echo "Ejecutando Híbrido con P=$P, T=$T (Total: $TOTAL_CORES)..."
            TIME=$(OMP_NUM_THREADS=$T mpirun --oversubscribe -np $P ./hybrid_exec $M $N $T | grep "Total time" | awk '{print $4}')
            echo "hybrid,$P,$T,$TOTAL_CORES,$GRID_SIZE,$TIME" >> "$OUTPUT_FILE"
        fi
    done
done
echo ""

# --- Fase 3: Finalización ---
echo "--- ✅ Todas las pruebas canónicas han finalizado. ---"
echo "Puedes revisar los resultados en el archivo: $OUTPUT_FILE"
