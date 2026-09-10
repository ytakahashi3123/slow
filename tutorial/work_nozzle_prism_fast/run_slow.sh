#!/bin/bash

# slow をインストールしていない場合でも動くよう src/ を PYTHONPATH に加える
# (pip install -e . 済みなら `slow > $LOG` だけでよい)
SRC=$(cd "$(dirname "$0")/../../src" && pwd)
LOG=log_slow

time PYTHONPATH="$SRC:$PYTHONPATH" python3 -m slow > $LOG

# numba が効いているかはログの先頭に出る
grep "numerical kernels" $LOG
