#!/usr/bin/env bash

TEX=(pdflatex -interaction=batchmode -shell-escape)
# BIB=(bibtex -terse)
BIB=(biber --nolog -m 99)
# -m 99 suppresses cross references that are not directly cited

DEP=(extra.tex)

NAME=''
if [ -z "$NAME" ]; then
    NAME=()
    for tex in *.tex; do
        for dep in "${DEP[@]}"; do
            [ "$tex" == "$dep" ] && continue 2
        done
        NAME+=("$tex")
    done

    if [ ${#NAME[@]} -ne 1 ]; then
        echo "Found ${#NAME[@]} .tex files"
        for x in "${NAME[@]}"; do
            echo "$x"
        done
        exit 1
    fi
    NAME="${NAME[0]%.tex}"
fi

if [ "$1" == "clean" ]; then
    clean=()
    for ext in pdf aux log out toc lof lot bbl bcf blg run.xml nav snm; do
        clean+=("${NAME}.${ext}")
    done
    [ ${#clean[@]} -gt 0 ] && rm -fv "${clean[@]}"
    exit
fi

PDF="${NAME}.pdf"

if [ "$1" == "optimize" ]; then
    if [ -f "${PDF}" ]; then
        for gs in gs gswin64; do
            if command -v "$gs" 2>&1 >/dev/null; then
                "$gs" -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 \
                    -dNOPAUSE -dQUIET -dBATCH -sOutputFile="tmp_${PDF}" "${PDF}"
                mv "tmp_${PDF}" "${PDF}"
                break
            fi
        done
    fi
    exit
fi

DEP=("${NAME}.tex" "${DEP[@]}")

md5() { md5sum "${NAME}.$1" 2> /dev/null; }

newerDEPS() {
    for x in "${DEP[@]}"; do
        [ "$x" -nt "$PDF" ] && return 0
    done
    return 1
}

warn=0 # show warnings: 0 = on error, 1 = always
for (( i=1, n=1; i<=n; ++i )); do
    md5_aux="$(md5 aux)"
    md5_bcf="$(md5 bcf)"
    # run LaTeX
    if (( i != 1 )) || newerDEPS; then
        printf "\e[32;1m$i\e[0m\n"
        if ! "${TEX[@]}" "${NAME}" > /dev/null; then
            warn=1
            break
        fi
    fi
    # check if need to run multiple times
    if (( i == 1 )) && ( # update bibliography
        newerDEPS || [ "$md5_bcf" != "$(md5 bcf)" ]
    ); then
        printf '\e[32;1mbib\e[0m\n'
        "${BIB[@]}" "${NAME}" | awk '''
            sub(/^WARN/,"\033[33m&\033[0m") || \
            sub(/^ERROR/,"\033[31m&\033[0m") \
            { print }
        '''
        [ ${PIPESTATUS[0]} -eq 0 ] || break
        (( ++n ))
    elif [ "$md5_aux" != "$(md5 aux)" ]; then # if aux file updated
        (( ++n ))
    fi
done
if (( warn != 0 )); then
    awk '''
        sub(/.*Warning:/,"\033[33m&\033[0m") || \
        sub(/^!.*/,"\033[31m&\033[0m") { p=1 }
        /^$/ { p=0 }
        p # print if p != 0
    ''' "${NAME}.log"
fi
