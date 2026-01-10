#!/usr/bin/env bash

set -e

cd "${0%/*}"
mkdir -p build

IFS=',' read -ra compilers <<< "$1"
[ ${#compilers[@]} -eq 0 ] && compilers=(g++ clang++)

IFS=',' read -ra stds <<< "$2"
[ ${#stds[@]} -eq 0 ] && stds=(c++20 c++23)

IFS=',' read -ra suite <<< "$3"
[ ${#suite[@]} -eq 0 ] && suite=(tests*.cpp)

srcs=(test.hpp ../include/*.hpp)
exes=()
pids=()

for tests in "${suite[@]}"; do
    newest="test.sh"
    for src in "$tests" "${srcs[@]}"; do
        if [ "$src" -nt "$newest" ]; then
            newest="$src"
        fi
    done

    for comp in "${compilers[@]}"; do
    for std in "${stds[@]}"; do
        exe="build/${tests%%.cpp}-$std-$comp"
        exes+=("$exe")
        if [ ! -f "$exe" ] || [ "$newest" -nt "$exe" ]; then
            if [ "$comp" == 'cl' ]; then
                cmd=("$comp" /EHsc /std:"$std" /O2 /W4 /WX /wd4702 \
                    /I. /I../include \
                    /Fe"$exe" "${tests}")
            else
                cmd=("$comp" -std="$std" -O3 -Wall -Wextra -Werror -pedantic \
                    -I. -I../include \
                    "${tests}" -o "$exe")
            fi
            echo -e '\033[34m'"${cmd[@]}"'\033[0m'
            if [ "$comp" == 'cl' ]; then
                "${cmd[@]}"
            else
                "${cmd[@]}" &
            fi
            pids+=($!)
        fi
    done
    done
done

for pid in ${pids[@]}; do
    wait $pid || exit 1
done

for exe in "${exes[@]}"; do
    echo "$exe"
    ./"$exe"
done
