#!/bin/bash
# Download and extract all remaining .d.rar files for PXD069059
# Already have: 22-Prot-1849_BD5_1_11838, 22-Prot-1849_BD5_1_11876, 22-Prot-1852_BD8_1_11874

set -e
DEST="/globalscratch/dateschn/claudius-data/raw/PXD069059"
UNRAR="/home/dateschn/unrar/unrar"
FTP_BASE="ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2026/01/PXD069059"

# All 67 .d.rar files
FILES=(
22-Prot-1849_BD5_1_11838.d.rar
22-Prot-1849_BD5_1_11876.d.rar
22-Prot-1850_BD6_1_11839.d.rar
22-Prot-1850_BD6_1_11872.d.rar
22-Prot-1851_BD7_1_11840.d.rar
22-Prot-1851_BD7_1_11903.d.rar
22-Prot-1852_BD8_1_11841.d.rar
22-Prot-1852_BD8_1_11874.d.rar
22-Prot-1853_BE1_1_11843.d.rar
22-Prot-1853_BE1_1_11878.d.rar
22-Prot-1854_BE2_1_11844.d.rar
22-Prot-1854_BE2_1_11879.d.rar
22-Prot-1855_BE3_1_11845.d.rar
22-Prot-1855_BE3_1_11880.d.rar
22-Prot-1856_BE4_1_11846.d.rar
22-Prot-1856_BE4_1_11881.d.rar
22-Prot-1857_BE5_1_11848.d.rar
22-Prot-1857_BE5_1_11883.d.rar
22-Prot-1858_BE6_1_11849.d.rar
22-Prot-1858_BE6_1_11884.d.rar
22-Prot-1859_BE7_1_11850.d.rar
22-Prot-1859_BE7_1_11885.d.rar
22-Prot-1860_BE8_1_11851.d.rar
22-Prot-1860_BE8_1_11886.d.rar
22-Prot-1861_GA1_1_11853.d.rar
22-Prot-1861_GA1_1_11888.d.rar
22-Prot-1862_GA2_1_11854.d.rar
22-Prot-1862_GA2_1_11889.d.rar
22-Prot-1863_GA3_1_11855.d.rar
22-Prot-1863_GA3_1_11890.d.rar
22-Prot-1864_GA4_1_11856.d.rar
22-Prot-1864_GA4_1_11891.d.rar
22-Prot-1865_GA5_1_11858.d.rar
22-Prot-1865_GA5_1_11893.d.rar
22-Prot-1866_GA6_1_11862.d.rar
22-Prot-1866_GA6_1_11894.d.rar
22-Prot-1867_GA7_1_11863.d.rar
22-Prot-1867_GA7_1_11895.d.rar
22-Prot-1868_GA8_1_11864.d.rar
22-Prot-1868_GA8_1_11896.d.rar
22-Prot-1869_GB1_1_11866.d.rar
22-Prot-1869_GB1_1_11898.d.rar
22-Prot-1870_GB2_1_11867.d.rar
22-Prot-1870_GB2_1_11899.d.rar
22-Prot-1871_GB3_1_11868.d.rar
22-Prot-1871_GB3_1_11900.d.rar
22-Prot-1872_GB4_1_11869.d.rar
22-Prot-1872_GB4_1_11901.d.rar
23-Prot-1759_RA1_1_15700.d.rar
23-Prot-1760_RA2_1_15701.d.rar
23-Prot-1761_RA3_1_15703.d.rar
23-Prot-1762_RA4_1_15704.d.rar
23-Prot-1763_RA5_1_15706.d.rar
23-Prot-1764_RA6_1_15707.d.rar
23-Prot-1765_RA7_1_15709.d.rar
23-Prot-1766_RA8_1_15710.d.rar
23-Prot-1767_RB1_1_15712.d.rar
23-Prot-1768_RB2_1_15713.d.rar
23-Prot-1769_RB3_1_15715.d.rar
23-Prot-1770_RB4_1_15716.d.rar
23-Prot-1771_RB5_1_15718.d.rar
23-Prot-1772_RB6_1_15719.d.rar
23-Prot-1773_RB7_1_15722.d.rar
23-Prot-1774_RB8_1_15723.d.rar
23-Prot-1775_RC1_1_15728.d.rar
23-Prot-1776_15ul_RC2_1_15729.d.rar
23-Prot-1777_RC3_1_15725.d.rar
)

DOWNLOADED=0
EXTRACTED=0
SKIPPED=0
FAILED=0
TOTAL=${#FILES[@]}

echo "=== PXD069059 Full Download ==="
echo "Total files: $TOTAL"
echo "Destination: $DEST"
echo "Started: $(date)"
echo ""

for f in "${FILES[@]}"; do
    d_name="${f%.rar}"  # e.g. 22-Prot-1849_BD5_1_11838.d
    d_dir="$DEST/$d_name"
    rar_file="$DEST/$f"

    # Skip if .d directory already exists and has analysis.tdf
    if [ -f "$d_dir/analysis.tdf" ]; then
        SKIPPED=$((SKIPPED + 1))
        echo "[$((SKIPPED + DOWNLOADED + FAILED))/$TOTAL] SKIP $d_name (already extracted)"
        continue
    fi

    # Download if .rar doesn't exist
    if [ ! -f "$rar_file" ]; then
        echo "[$((SKIPPED + DOWNLOADED + FAILED + 1))/$TOTAL] Downloading $f ..."
        if ! wget -q --timeout=300 --tries=3 "$FTP_BASE/$f" -O "$rar_file"; then
            echo "  FAILED to download $f"
            FAILED=$((FAILED + 1))
            rm -f "$rar_file"
            continue
        fi
        DOWNLOADED=$((DOWNLOADED + 1))
        echo "  Downloaded ($(du -h "$rar_file" | cut -f1))"
    else
        echo "[$((SKIPPED + DOWNLOADED + FAILED + 1))/$TOTAL] RAR exists: $f"
    fi

    # Extract
    echo "  Extracting $f ..."

    # Check if archive has .d parent directory
    has_parent=true
    while IFS= read -r entry; do
        if [ -n "$entry" ] && [ "$entry" != "$d_name" ] && [[ ! "$entry" == "$d_name/"* ]]; then
            has_parent=false
            break
        fi
    done < <("$UNRAR" lb "$rar_file" 2>/dev/null)

    if [ "$has_parent" = true ]; then
        "$UNRAR" x -o+ "$rar_file" "$DEST/" > /dev/null 2>&1
    else
        mkdir -p "$d_dir"
        "$UNRAR" x -o+ "$rar_file" "$d_dir/" > /dev/null 2>&1
    fi

    if [ -f "$d_dir/analysis.tdf" ]; then
        EXTRACTED=$((EXTRACTED + 1))
        echo "  Extracted OK ($(du -sh "$d_dir" | cut -f1))"
    else
        echo "  WARNING: extraction may have failed for $d_name"
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "=== Summary ==="
echo "Downloaded: $DOWNLOADED"
echo "Extracted: $EXTRACTED"
echo "Skipped (already done): $SKIPPED"
echo "Failed: $FAILED"
echo "Finished: $(date)"

# Count total .d directories
N_DIRS=$(find "$DEST" -maxdepth 1 -name "*.d" -type d | wc -l)
echo "Total .d directories: $N_DIRS / $TOTAL"
