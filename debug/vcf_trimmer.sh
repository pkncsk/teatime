#!/bin/bash
# Usage: ./trim_gnomad_info.sh input.vcf.gz region output.vcf.gz

input_vcf="$1"
region="$2"
output_vcf="$3"

tabix "$input_vcf" "$region" | \
awk 'BEGIN{OFS="\t"}
    /^##INFO=/ {
        if ($0 ~ /ID=AC/ || $0 ~ /ID=AN/ || $0 ~ /ID=AF/)
            print
        next
    }
    /^##/ { print; next }
    /^#/ {
        print
        next
    }
    {
        split($8, info_fields, ";")
        new_info = ""
        for(i in info_fields){
            if (info_fields[i] ~ /^AC=/ || info_fields[i] ~ /^AN=/ || info_fields[i] ~ /^AF=/){
                if (new_info != "") new_info = new_info ";"
                new_info = new_info info_fields[i]
            }
        }
        $8 = new_info
        print
    }
' | bgzip > "$output_vcf"

tabix -p vcf "$output_vcf"
