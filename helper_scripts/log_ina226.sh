#!/usr/bin/env bash
# log_ina226.sh  stream INA226 rails to CSV until Ctrl-C
# Usage: sudo ./log_ina226.sh [output.csv] [interval_sec]
#   output.csv    (default: ./ina226_YYYYmmdd_HHMMSS.csv)
#   interval_sec  (default: 0.1)
#
# THIS SCRIPT IS MEANT TO BE RUN ON THE VERSAL!
#
# Examples:
# ./log_ina226.sh
# ./log_ina226.sh rails.csv 0.5

set -u # (no -e)

outfile="${1:-"ina226_$(date +%Y%m%d_%H%M%S).csv"}"
interval="${2:-0.1}"

# --- discover INA226 hwmon devices ---
declare -a HWMONS=()
for d in /sys/class/hwmon/hwmon*; do
	[[ -f "$d/name" ]] || continue
	if grep -q '^ina226' "$d/name"; then HWMONS+=("$d"); fi
done
if [[ ${#HWMONS[@]} -eq 0 ]]; then
	echo "No INA226 hwmon devices found under /sys/class/hwmon" >&2
	exit 1
fi

# Build columns: LABEL_mV, LABEL_mA, LABEL_uW
declare -a LABELS=() MV=() MA=() UW=()
for d in "${HWMONS[@]}"; do
	lbl=""
	[[ -f "$d/label" ]] && lbl="$(tr -d '\n' <"$d/label")"
	if [[ -z "$lbl" ]]; then
		devlink="$(readlink -f "$d/device" || true)" # e.g. .../i2c-5/5-0041
		lbl="$(basename "$devlink" 2>/dev/null || echo hwmon)"
	fi
	lbl="${lbl// /_}"
	lbl="${lbl//[^A-Za-z0-9_]/_}"

	mvf="$d/in1_input"    # bus voltage (mV)
	maf="$d/curr1_input"  # current (mA)
	uwf="$d/power1_input" # power (µW)
	if [[ -r "$mvf" && -r "$maf" && -r "$uwf" ]]; then
		LABELS+=("$lbl")
		MV+=("$mvf")
		MA+=("$maf")
		UW+=("$uwf")
	fi
done
if [[ ${#LABELS[@]} -eq 0 ]]; then
	echo "INA226 devices found, but no readable in1/curr1/power1 files." >&2
	exit 1
fi

# Respect sensor update_interval (max across devices)
max_upd_ms=0
for d in "${HWMONS[@]}"; do
	if [[ -r "$d/update_interval" ]]; then
		ui="$(cat "$d/update_interval" 2>/dev/null || echo 0)"
		[[ "$ui" =~ ^[0-9]+$ ]] || ui=0
		((ui > max_upd_ms)) && max_upd_ms="$ui"
	fi
done

# Compute sleep in microseconds, using the larger of requested vs update_interval
req_us=$(awk -v s="$interval" 'BEGIN{ printf("%d", s*1000000) }')
upd_us=$((max_upd_ms * 1000))
((req_us < upd_us)) && sleep_us="$upd_us" || sleep_us="$req_us"

have_usleep=0
command -v usleep >/dev/null 2>&1 && have_usleep=1
sleep_int=$(((sleep_us + 500000) / 1000000)) # round to nearest sec for fallback

# Write header
{
	printf "ts_ms"
	for ((i = 0; i < ${#LABELS[@]}; i++)); do
		printf ",%s_mV,%s_mA,%s_uW" "${LABELS[i]}" "${LABELS[i]}" "${LABELS[i]}"
	done
	printf "\n"
} >"$outfile"

rows=0
now_ms() { date +%s%3N 2>/dev/null || printf '%s000' "$(date +%s)"; }

echo "Logging to: $outfile  (interval: ${sleep_us}us)  Press Ctrl-C to stop."
trap 'echo; echo "Stopped. Rows written: $rows"; exit 0' INT TERM

# --- infinite loop ---
while :; do
	ts="$(now_ms)"
	line="$ts"
	for ((i = 0; i < ${#LABELS[@]}; i++)); do
		v_mv="$(cat "${MV[i]}" 2>/dev/null || echo)"
		v_ma="$(cat "${MA[i]}" 2>/dev/null || echo)"
		v_uw="$(cat "${UW[i]}" 2>/dev/null || echo)"
		line+=",$v_mv,$v_ma,$v_uw"
	done
	echo "$line" >>"$outfile"
	((rows++))

	if ((have_usleep)); then
		usleep "$sleep_us" || sleep "$sleep_int"
	else
		sleep "$sleep_int"
	fi
done
