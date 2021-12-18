#!/bin/bash
module load samtools
module load parallel

export refpath="/blue/jones/anujsharma/reference/Capsicum_annuum/Zunla1"
export refline="Zunla1"
export filename="$1"
export path="/blue/jones/anujsharma/scripts"
export gene="$2"
export filter="$3"

export SIZE=10
export SIZE2=75

REformat() {

	subnwrite() {
		pat=`awk '{ gsub("W","[AT]"); gsub("S","[CG]"); gsub("M","[AC]"); gsub("K","[GT]"); gsub("R","[AG]"); gsub("Y","[CT]"); gsub("B","[CGT]"); gsub("D","[AGT]"); gsub("H","[ACT]"); gsub("V","[ACG]"); gsub("N","."); gsub("Z",""); print }' <<< "$2"`
		printf "%b\t%b\t%b\n" "$1" "$pat" "$3" >> $path/REdata_formatted
	}
	
	while IFS="()/," read -ra chars
		do 
			if [ ${#chars[@]} -eq 2 ]
			then
				subnwrite $1 ${chars[0]}${chars[1]} ${#chars[0]}
			else
				[[ ${#chars[@]} -eq 3 ]] && pat1=${chars[0]} && pos1=$((chars[1]+${#chars[0]})) && pos2=$((-chars[2]))
				[[ ${#chars[@]} -eq 5 ]] && pat1=${chars[2]} && pos1=`echo {((-chars[0]))},${((chars[3]+${#chars[0]}))}` && pos2=`echo ${((chars[1]+${#chars[0]}))},${((-chars[4]))}`
				pat2=`tr "[ATGCUWSMKRYBDHV]" "]TACGAWSKMYRVHDB[" <<< $pat1 | rev`
				subnwrite $1 $pat1 $pos1
				subnwrite $1 $pat2 $pos2
			fi
		done <<< "$2"
}

[[ ! -f $path/REdata_formatted ]] && echo "Formatting database..." && while IFS=$'\t' read -r NAME PAT; do REformat $NAME ${PAT^^} ; done < $path/REdatabase

filtertarget() {
	
	IFS=':' read -r chr_targ pos_targ <<<"$filter"
	IFS='-' read -r pos_start pos_end <<< "$pos_targ"
	[[ -z $pos_end ]] && pos_start=0 && pos_end=9999999999
	
	TMPDIR=$(mktemp -d)
	awk 'BEGIN{OFS="\t"} { if ( $1 ~ /chr/ || ($1=="'$chr_targ'" && $2>'$pos_start' && $2<'$pos_end' )) print $0}' $filename > $TMPDIR/targlist.txt
	filename=$TMPDIR/targlist.txt
}

export -f filtertarget
[[ ! -z $filter ]] && filtertarget

findsite() {
	
	IFS=$'\t' read -r -a LINE <<< "$1"
	num=${#LINE[@]}
	CHR=${LINE[0]}
	POS=${LINE[1]}
	REF=${LINE[2]}
	ALT=${LINE[3]^^}
	MUT=${LINE[4]}
	B1R=${LINE[$((num-4))]}
	B1A=${LINE[$((num-3))]}
	B2R=${LINE[$((num-2))]}
	B2A=${LINE[$((num-1))]}
	
	[[ ${#REF} -gt 40 || ${#ALT} -gt 40 ]] && COMMENT='Possible_PCR_marker' || COMMENT=''

	OLD=`samtools faidx $refpath/$refline.chrs.fa $CHR:$(($POS-$SIZE2))-$(($POS+${#REF}-1+$SIZE2)) | grep -v '>' | tr -d '[:space:]' | tr a-z A-Z`
	NEW=${OLD:0:$SIZE2}$ALT${OLD:$SIZE2+${#REF}:$SIZE2}
	SAMF=${OLD:65:10}"["$REF"/"$ALT"]"${OLD:$SIZE2+${#REF}:-65}
	
	REO=()
	REN=()
	
	while IFS=$'\t' read -r NAME PAT CUT
		do
			[[ `grep -Eo "$PAT" <<< "$OLD"` && ! `grep "$PAT" <<< "$NEW"` ]] && REO+=($NAME)
			[[ ! `grep -Eo "$PAT" <<< "$OLD"` && `grep "$PAT" <<< "$NEW"` ]] && REN+=($NAME)
		done < $path/REdata_formatted_6mer ###

	OLD=${OLD:65:-65}
	NEW=${NEW:65:-65}
	
	printf "%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\n" "$CHR" "$POS" "$REF" "$ALT" "$MUT" "$B1R" "$B1A" "$B2R" "$B2A" "$OLD" "$NEW" "$SAMF" $(IFS=','; echo "${REO[*]}"$"\t""${REN[*]/%/*}") $COMMENT >> $gene.markers.txt ###
}

printf "%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\t%b\n" "chr" "pos"  "ref" "alt" "mutation" "R.ref" "R.alt" "S.ref" "S.alt" "ref.seq" "alt.seq" "samf" "ref.RE" "alt.RE" "misc" > $2.markers.txt ###

export -f findsite 
[[ `grep 'chr' <<< head -n1 $1` || `grep 'chr' <<< head -n1 $1` ]] && header=2 || header=1
tail -n +$header $filename | parallel --will-cite findsite ###

sort -k2 -n $gene.markers.txt > $gene.markers.sort.txt

[[ ! -z TMPDIR ]] && rm -rf $TMPDIR
# sbatch --wrap "/blue/jones/anujsharma/scripts/markers_wip gene.seg_SNPs.txt gene chr:start-end" --account "plantpath" --qos "plantpath" --time "01:00:00" --job-name "gene.marker" --mail-user "anujsharma@ufl.edu" --ntasks 1 --cpus-per-task 2 --mem 8g #Don't forget to change refline
