#!/bin/bash

# ==============================================================================
# Script Final para Generar Datos de Rendimiento para Gráficas de Escalabilidad
# ==============================================================================
#
# Este script automatiza la compilación, ejecución y análisis de rendimiento
# para tres versiones de un programa (MPI, OpenMP, Híbrido).
#

# Salir inmediatamente si un comando falla.
set -e

# --- CONFIGURACIÓN ---
# Rutas a los archivos de código fuente
MPI_SRC="Distribuido/marching_squares.cpp"
OMP_SRC="Compartido/marching_squares.cpp" # Usamos el de 'Compartido' para OMP
HYBRID_SRC="Hibrido/marching_squares.cpp"
SEQUENTIAL_SRC="Secuencial/marching_square.cpp" # El secuencial puro

# Nombres de los ejecutables de salida
MPI_EXE="Distribuido/marching_squares_mpi"
OMP_EXE="Compartido/marching_squares_omp"
HYBRID_EXE="Hibrido/marching_squares_hybrid"
SEQUENTIAL_EXE="Secuencial/marching_squares_sequential"

# Compiladores
MPI_COMPILER="mpicxx"
CXX_COMPILER="g++-14"

# Banderas de compilación
INCLUDE_PATH="/opt/homebrew/include"
LIBRARY_PATH="/opt/homebrew/lib"
CXX_FLAGS="-O3 -I${INCLUDE_PATH} -L${LIBRARY_PATH}"
OMP_FLAGS="-fopenmp"

# Parámetros de la simulación
GRID_SIZE_M=10000
GRID_SIZE_N=10000
LOG_FILE="performance_results.log"
# --- FIN DE LA CONFIGURACIÓN ---

print_header() {
    echo ""
    echo "##############################################################################"
    echo "## $1"
    echo "##############################################################################"
    echo ""
}

# Función robusta para extraer el tiempo de una salida estandarizada.
extract_time_ms() {
    grep 'FINAL_TIME_MS:' | cut -d':' -f2 | tr -d ' '
}

compile_all() {
    print_header "PASO 1: COMPILANDO TODOS LOS PROGRAMAS"
    echo "Compilando Secuencial..."
    ${CXX_COMPILER} ${CXX_FLAGS} ${SEQUENTIAL_SRC} -o ${SEQUENTIAL_EXE}
    echo "Compilando MPI puro..."
    ${MPI_COMPILER} ${CXX_FLAGS} ${MPI_SRC} -o ${MPI_EXE}
    echo "Compilando OpenMP puro..."
    ${MPI_COMPILER} ${CXX_FLAGS} ${OMP_FLAGS} ${OMP_SRC} -o ${OMP_EXE}
    echo "Compilando Híbrido..."
    ${MPI_COMPILER} ${CXX_FLAGS} ${OMP_FLAGS} ${HYBRID_SRC} -o ${HYBRID_EXE}
    echo "Compilación completada."
}

get_sequential_baseline() {
    print_header "PASO 2: OBTENIENDO TIEMPO BASE SECUENCIAL T(1)"
    local time_output=$(./${SEQUENTIAL_EXE} ${GRID_SIZE_M} ${GRID_SIZE_N})
    echo "${time_output}"
    local baseline_time=$(echo "${time_output}" | extract_time_ms)
    
    if [ -z "${baseline_time}" ]; then
        echo "ERROR: No se pudo extraer el tiempo base T(1)." >&2
        exit 1
    fi
    
    echo ""
    echo "TIEMPO BASE T(1) = ${baseline_time} ms"
    echo "${baseline_time}"
}

run_performance_tests() {
    local test_name="$1"
    local command_template="$2"
    local baseline_time="$3"
    local cores_to_test=(1 2 4 8)

    print_header "PASO: EJECUTANDO PRUEBAS DE ESCALABILIDAD PARA: ${test_name}"
    printf "%-10s, %-15s, %-15s, %-15s\n" "Cores(p)" "Tiempo T(p)[ms]" "Speedup S(p)" "Eficiencia E(p)[%]"

    for p in "${cores_to_test[@]}"; do
        local command_to_run="${command_template//\{p\}/$p}"
        
        local output
        if [[ "$test_name" == *"MPI"* || "$test_name" == *"HÍBRIDO"* ]] && [ "$p" -eq 1 ]; then
            command_to_run="${command_to_run//mpirun -np 1 /}"
        fi
        
        output=$(eval "${command_to_run}")
        local time_p=$(echo "${output}" | extract_time_ms)

        if [ -z "${time_p}" ] || (( $(echo "$time_p <= 0" | bc -l) )); then
            printf "%-10d, %-15s, %-15s, %-15s\n" "$p" "ERROR" "ERROR" "ERROR"
            continue
        fi

        local speedup=$(echo "scale=4; ${baseline_time} / ${time_p}" | bc)
        local efficiency=$(echo "scale=4; (${speedup} / ${p}) * 100" | bc)
        printf "%-10d, %-15.3f, %-15.3f, %-15.2f\n" "$p" "$time_p" "$speedup" "$efficiency"
    done
}

main() {
    rm -f "${LOG_FILE}"
    exec &> >(tee -a "${LOG_FILE}")
    echo "Iniciando script de pruebas de rendimiento - $(date)"
    
    compile_all
    local T1=$(get_sequential_baseline)

    local mpi_cmd_template="mpirun -np {p} ./${MPI_EXE} ${GRID_SIZE_M} ${GRID_SIZE_N}"
    local omp_cmd_template="./${OMP_EXE} ${GRID_SIZE_M} ${GRID_SIZE_N} {p}"
    local hybrid_cmd_template="mpirun -np {p} ./${HYBRID_EXE} ${GRID_SIZE_M} ${GRID_SIZE_N} 1"

    run_performance_tests "MPI PURO" "${mpi_cmd_template}" "${T1}"
    run_performance_tests "OPENMP PURO" "${omp_cmd_template}" "${T1}"
    run_performance_tests "HÍBRIDO (1 Hilo/Proceso)" "${hybrid_cmd_template}" "${T1}"

    print_header "DATOS PARA GRÁFICAS GENERADOS"
    echo "Copia las tablas de arriba a tu programa de hojas de cálculo."
}

main

