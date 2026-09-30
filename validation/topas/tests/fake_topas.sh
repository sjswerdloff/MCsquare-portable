#!/bin/bash
# Stand-in for the topas binary in the wrapper tests. No transport: writes outputs shaped by FAKE_MODE.
# Reads the grid from ./stage1_base.txt and the EM module from ./run.txt, as the real run would.
set -u
mode=${FAKE_MODE:-ok}
bins() { awk -F= -v k="i:Ge/Phantom/$1" '{g=$1; gsub(/[ \t]/,"",g)} g==k {v=$2; gsub(/[ \t]/,"",v); print v}' stage1_base.txt; }
n=$(( $(bins XBins) * $(bins YBins) * $(bins ZBins) * 8 ))
icru=0; grep -q '"g4em-standard_opt4"' run.txt && icru=1
[ "$mode" = wrongem ] && icru=$(( 1 - icru ))
[ "$mode" = slow ] && sleep 2
[ "$mode" = noem ] || echo "Use ICRU90 data                                    $icru"
write() { head -c "$2" /dev/zero > "$1"; }
case "$mode" in
    exit3) exit 3 ;;
esac
write dose.bin "$n"
[ "$mode" = missing ] || write dose_all.bin "$n"
[ "$mode" = wrongsize ] && write dose_all.bin $(( n - 8 ))
echo "# header" > dose.binheader
if [ "$mode" = empty ]; then : > dose_all.binheader; else echo "# header" > dose_all.binheader; fi
exit 0
