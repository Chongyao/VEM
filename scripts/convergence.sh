#!/bin/bash

# 设置路径
BUILD_DIR="./build"
BIN_DIR="$BUILD_DIR/bin"
SCRIPT_DIR="./scripts"
OUTPUT_DIR="./results"
# 确保编译是新的
echo "--- Compiling Project ---"
cmake -B $BUILD_DIR -S . -DCMAKE_BUILD_TYPE=Release
cmake --build $BUILD_DIR --target test_beam_gravity -j4

# 创建结果目录
mkdir -p $OUTPUT_DIR

# 定义网格密度列表
SEEDS=(100 200 500 1000 2000)
LOG_FILE="$OUTPUT_DIR/convergence_log.txt"
rm -f $LOG_FILE

echo "--- Starting Convergence Test ---"
echo "Num_Elements L2_Error" > $LOG_FILE

for N in "${SEEDS[@]}"; do
    MESH_FILE="$OUTPUT_DIR/beam_${N}.vtu"
    
    echo "Processing N=$N ..."
    
    # 1. 生成网格
    python3 $SCRIPT_DIR/cvt.py --n_seeds $N --output $MESH_FILE > /dev/null
    
    # 2. 运行求解器并提取结果
    # 我们grep抓取 "[RESULT]" 开头的行
    OUTPUT=$($BIN_DIR/test_beam_gravity $MESH_FILE)
    RESULT_LINE=$(echo "$OUTPUT" | grep "\[RESULT\]")
    
    if [ -z "$RESULT_LINE" ]; then
        echo "  Error: Simulation failed for N=$N"
        echo "$OUTPUT" # 打印错误日志
    else
        # 提取数字部分: 200 0.012345
        DATA=$(echo "$RESULT_LINE" | awk '{print $2, $3}')
        echo "  Done. Result: $DATA"
        echo "$DATA" >> $LOG_FILE
    fi
done

echo "--- Convergence Test Finished ---"
cat $LOG_FILE

# (可选) 调用 Python 画图
# python3 $SCRIPT_DIR/plot_error.py $LOG_FILE
