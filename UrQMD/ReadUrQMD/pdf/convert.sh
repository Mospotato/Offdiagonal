#!/bin/bash
for input in *.pdf; do
    # Extract the base filename (without extension)
    output="${input%.*}.png"

    # Convert PDF to PNG using pdftoppm
    pdftoppm -png -r 300 -scale-to 0 -cropbox -singlefile "$input" "${input%.*}"

    # OR Uncomment to use ImageMagick
    # convert -density 300 "$input" "$output"

    # OR Uncomment to use Ghostscript
    # gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r300 -sOutputFile="$output" "$input"

    echo "Converted $input to $output"
done
