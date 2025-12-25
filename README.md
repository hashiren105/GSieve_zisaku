### GSieve_zisaku

ビルド方法：
```
make clean && make
```

実行方法：
```
./obj/gsieve "INPUT_FILE"
```

うまく動かないときはシグマなどの設定を調整してください。

### GH列挙スクリプト

1. `python3.10 -m pip install --user fpylll` で `fpylll` を入れます。
2. 次のようにスクリプトを走らせると、GH（ガウスヒューリスティック）比が近いベクトルを探します。

```
python3.10 scripts/enum_fpylll.py input/ideallatticeindex23seed0.txt \
    --factor 1.05 --max-factor 1.2 --step 0.01
```

`--factor`〜`--max-factor` の範囲で試行し、目的の比率に届くとそのベクトルを表示して終了します。必要に応じて `--step` や `--show`（出力件数）を調整してください。
