#!/usr/bin/env bash
#
# cortado batch runner: run measure_pahs_general.py (cortado) on every source
# listed in a table.
#
# Usage:
#   run_pah_fits.sh [options] TABLE [-- extra measure_pahs_general.py arguments]
#
# Options:
#   -n, --dry-run        print the commands without running them
#   -o, --outdir DIR     output directory for all fits (default: ./pah_fits_table)
#   -p, --python EXE     Python interpreter (default: $PYTHON, else python)
#   -f, --fitter PATH    path to measure_pahs_general.py (default: next to this script)
#   -m, --model-range R  model drawn in the figures: 'full' (default) from the bluest
#                        data to the red edge of the reddest filter, dashed outside
#                        the spectra; 'data' only over the spectra
#   -h, --help           show this help
#
# Anything after "--" is passed to every fit, e.g.
#   run_pah_fits.sh sources.csv -- --partial-pahs --flux-unit Jy
#
# TABLE format: comma-separated, one header row. Blank lines and lines starting
# with '#' are ignored. Fields are split on commas only (no quoting), so values
# must not contain commas; spaces at either end of a field are removed.
#
#   id      source name, used for the output files              (required)
#   spec    LRS spectrum, or a combined LRS+IRS file             (required)
#   z       redshift                                             (required)
#   irs     separate IRS spectrum                                (optional)
#   seam    split a combined file: lam_obs > seam (um) is IRS    (optional)
#   args    extra fitter arguments for this source only,         (optional)
#           space-separated, e.g. "--continuum poly"
#
# Every other column BAND that has a partner column BAND_err is photometry,
# passed as --phot BAND:value:error. BAND is a filter alias known to the fitter
# (MIPS24, MIPS70, PACS70, PACS100, F770W, ...), an SVO filter ID (e.g.
# JWST/MIRI.F1800W), or a number = observed wavelength in um, fit as a delta
# function (the model value at exactly that wavelength; never matched to a
# filter). Values are in the same flux unit as the spectra
# (mJy unless --flux-unit is given after "--"). A band is skipped for a source
# when its value or error is empty, not a number (e.g. nan), or the error is
# not positive. Column names are case-insensitive except for BAND names,
# which are passed through as written.
#
# Relative paths in spec / irs are taken relative to the table's directory.
# Each fit writes <outdir>/<id>_features.ecsv (fluxes, luminosities, EWs, S/N
# per feature) and <outdir>/<id>.png: data (LRS black, IRS grey), photometry at
# its calibration wavelength with the filter half-max width, and the model (red;
# solid on the spectra, dashed where extrapolated; red crosses = model through
# each filter) with its continuum (blue), PAH (magenta) and line (green)
# components and the attenuation (dotted, right axis). The fitter's printout
# goes to <outdir>/logs/<id>.log. A failed source is reported and the run
# continues; the exit status is the number of failures.

set -u

usage() { sed -n '3,/^$/s/^# \{0,1\}//p' "$0"; }

script_dir=$(cd "$(dirname "$0")" && pwd)
fitter="$script_dir/measure_pahs_general.py"
python_exe="${PYTHON:-python}"
outdir="pah_fits_table"
model_range="full"
dry_run=0
table=""
extra=()

need_value() { [ $# -ge 2 ] && [ -n "$2" ] || { echo "option $1 needs a value" >&2; exit 2; }; }

while [ $# -gt 0 ]; do
	case "$1" in
		-n|--dry-run)     dry_run=1 ;;
		-o|--outdir)      need_value "$@"; outdir="$2"; shift ;;
		-p|--python)      need_value "$@"; python_exe="$2"; shift ;;
		-f|--fitter)      need_value "$@"; fitter="$2"; shift ;;
		-m|--model-range) need_value "$@"; model_range="$2"; shift ;;
		-h|--help)        usage; exit 0 ;;
		--)           shift; extra=("$@"); break ;;
		-*)           echo "unknown option: $1" >&2; usage >&2; exit 2 ;;
		*)            if [ -z "$table" ]; then table="$1"; else echo "unexpected argument: $1" >&2; exit 2; fi ;;
	esac
	shift
done

[ -n "$table" ] || { usage >&2; exit 2; }
[ -f "$table" ] || { echo "table not found: $table" >&2; exit 2; }
[ -f "$fitter" ] || { echo "fitter not found: $fitter (use --fitter)" >&2; exit 2; }
case "$model_range" in full|data) ;; *) echo "--model-range must be 'full' or 'data', not '$model_range'" >&2; exit 2 ;; esac
table_dir=$(cd "$(dirname "$table")" && pwd)

trim() {   # strip surrounding spaces, tabs and a Windows carriage return
	local s="$1"
	s="${s%$'\r'}"
	s="${s#"${s%%[![:space:]]*}"}"
	s="${s%"${s##*[![:space:]]}"}"
	printf '%s' "$s"
}
is_number() { awk -v x="$1" 'BEGIN { exit !(x ~ /^[-+]?([0-9]+[.]?[0-9]*|[.][0-9]+)([eE][-+]?[0-9]+)?$/) }'; }
is_positive() { awk -v x="$1" 'BEGIN { exit !(x + 0 > 0) }'; }
lower() { printf '%s' "$1" | tr '[:upper:]' '[:lower:]'; }
resolve() {   # relative paths are relative to the table's directory
	case "$1" in /*|~*) printf '%s' "$1" ;; *) printf '%s/%s' "$table_dir" "$1" ;; esac
}
split_csv() {   # split $1 on commas into the global array 'fields', keeping empty fields
	local IFS=,
	read -r -a fields <<< "$1,__end__"
	unset 'fields[${#fields[@]}-1]'
}

# ---- header
header=""
while IFS= read -r line || [ -n "$line" ]; do
	line=$(trim "$line")
	[ -z "$line" ] && continue
	case "$line" in \#*) continue ;; esac
	header="$line"
	break
done < "$table"
[ -n "$header" ] || { echo "no header row in $table" >&2; exit 2; }

split_csv "$header"
names=()
for f in "${fields[@]}"; do names+=("$(trim "$f")"); done
ncol=${#names[@]}

i_id=-1; i_spec=-1; i_z=-1; i_irs=-1; i_seam=-1; i_args=-1
for ((k = 0; k < ncol; k++)); do
	case "$(lower "${names[k]}")" in
		id)   i_id=$k ;;
		spec) i_spec=$k ;;
		z)    i_z=$k ;;
		irs)  i_irs=$k ;;
		seam) i_seam=$k ;;
		args) i_args=$k ;;
	esac
done
for req in id spec z; do
	eval "idx=\$i_$req"
	[ "$idx" -ge 0 ] || { echo "table needs a '$req' column (header: $header)" >&2; exit 2; }
done

# photometry: columns BAND with a partner BAND_err
band_idx=(); err_idx=()
for ((k = 0; k < ncol; k++)); do
	case "$(lower "${names[k]}")" in id|spec|z|irs|seam|args|*_err) continue ;; esac
	for ((j = 0; j < ncol; j++)); do
		if [ "$(lower "${names[j]}")" = "$(lower "${names[k]}")_err" ]; then
			band_idx+=("$k"); err_idx+=("$j")
		fi
	done
done
for ((k = 0; k < ncol; k++)); do   # warn about columns that will be ignored
	case "$(lower "${names[k]}")" in id|spec|z|irs|seam|args) continue ;; esac
	used=0
	for j in ${band_idx[@]+"${band_idx[@]}"} ${err_idx[@]+"${err_idx[@]}"}; do [ "$j" -eq "$k" ] && used=1; done
	[ $used -eq 1 ] || echo "note: column '${names[k]}' ignored (photometry needs a pair BAND, BAND_err)" >&2
done
bands=""
for k in ${band_idx[@]+"${band_idx[@]}"}; do bands="$bands ${names[k]}"; done
echo "cortado: $("$python_exe" "$fitter" --version 2>/dev/null | sed 's/^cortado //')   table: $table   photometry bands:${bands:- none}"

# ---- rows
[ $dry_run -eq 1 ] || mkdir -p "$outdir/logs"
n_ok=0; n_fail=0; failed=""
seen_header=0
while IFS= read -r line || [ -n "$line" ]; do
	line=$(trim "$line")
	[ -z "$line" ] && continue
	case "$line" in \#*) continue ;; esac
	if [ $seen_header -eq 0 ]; then seen_header=1; continue; fi

	split_csv "$line"
	if [ ${#fields[@]} -ne $ncol ]; then
		echo "skipping row with ${#fields[@]} fields (header has $ncol): $line" >&2
		n_fail=$((n_fail + 1)); failed="$failed (bad row)"; continue
	fi
	vals=()
	for f in "${fields[@]}"; do vals+=("$(trim "$f")"); done

	id="${vals[i_id]}"; spec="${vals[i_spec]}"; z="${vals[i_z]}"
	if [ -z "$id" ] || [ -z "$spec" ] || ! is_number "$z"; then
		echo "skipping row with missing id/spec or non-numeric z: $line" >&2
		n_fail=$((n_fail + 1)); failed="$failed ${id:-(no id)}"; continue
	fi
	spec=$(resolve "$spec")
	if [ ! -f "$spec" ]; then
		echo "[$id] spectrum not found: $spec" >&2
		n_fail=$((n_fail + 1)); failed="$failed $id"; continue
	fi

	cmd=("$python_exe" "$fitter" "$spec" --z "$z" --name "$id" --outdir "$outdir" --model-range "$model_range")
	if [ $i_irs -ge 0 ] && [ -n "${vals[i_irs]}" ]; then
		irs=$(resolve "${vals[i_irs]}")
		[ -f "$irs" ] || { echo "[$id] IRS spectrum not found: $irs" >&2; n_fail=$((n_fail + 1)); failed="$failed $id"; continue; }
		cmd+=(--irs "$irs")
	fi
	if [ $i_seam -ge 0 ] && [ -n "${vals[i_seam]}" ]; then
		cmd+=(--seam "${vals[i_seam]}")
	fi
	if [ $i_args -ge 0 ] && [ -n "${vals[i_args]}" ]; then
		read -r -a row_args <<< "${vals[i_args]}"
		cmd+=("${row_args[@]}")
	fi
	n_phot=0
	for ((b = 0; b < ${#band_idx[@]}; b++)); do
		v="${vals[${band_idx[b]}]}"; e="${vals[${err_idx[b]}]}"
		if is_number "$v" && is_number "$e" && is_positive "$e"; then
			cmd+=(--phot "${names[${band_idx[b]}]}:$v:$e")
			n_phot=$((n_phot + 1))
		fi
	done
	cmd+=(${extra[@]+"${extra[@]}"})

	if [ $dry_run -eq 1 ]; then
		printf '%q ' "${cmd[@]}"; echo
		continue
	fi
	printf '[%s] z=%s, %d photometric bands ... ' "$id" "$z" "$n_phot"
	if "${cmd[@]}" > "$outdir/logs/$id.log" 2>&1; then
		echo "ok"
		n_ok=$((n_ok + 1))
	else
		echo "FAILED (see $outdir/logs/$id.log)"
		n_fail=$((n_fail + 1)); failed="$failed $id"
	fi
done < "$table"

if [ $dry_run -eq 0 ]; then
	echo "done: $n_ok fitted, $n_fail failed${failed:+ ->$failed}; outputs in $outdir"
fi
exit $n_fail
