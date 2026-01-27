"""
Rules for managing FASTA databases.

Handles:
- Linking local FASTA files
- Downloading from UniProt by proteome ID
- Adding contaminants (cRAP)
- Mapping datasets to organisms
"""

import os
from pathlib import Path


# Constrain wildcards to prevent ambiguous rule matching
wildcard_constraints:
    organism="[^/]+",  # organism cannot contain slashes
    accession="PXD[0-9]+"  # accession must be PRIDE format


def get_organism_for_dataset(accession):
    """Get organism name for a dataset from config."""
    metadata = config.get("dataset_metadata", {}).get(accession, {})
    return metadata.get("organism", "human")  # Default to human


def get_fasta_for_organism(organism):
    """Get FASTA path for an organism (local or to be downloaded)."""
    org_config = config.get("organisms", {}).get(organism, {})
    local = org_config.get("local_fasta")
    if local and Path(local).exists():
        return local
    return None


def organism_includes_contaminants(organism):
    """Check if organism's local FASTA already includes contaminants."""
    org_config = config.get("organisms", {}).get(organism, {})
    return org_config.get("includes_contaminants", False)


rule get_or_download_fasta:
    """
    Get FASTA for an organism - link local or download from UniProt.
    """
    output:
        fasta="resources/fasta/{organism}.fasta"
    params:
        local_fasta=lambda wildcards: config.get("organisms", {}).get(wildcards.organism, {}).get("local_fasta"),
        proteome_id=lambda wildcards: config.get("organisms", {}).get(wildcards.organism, {}).get("proteome_id"),
        taxon_id=lambda wildcards: config.get("organisms", {}).get(wildcards.organism, {}).get("taxon_id")
    log:
        "logs/fasta/{organism}.log"
    run:
        from pathlib import Path
        import urllib.request
        import shutil

        output_path = Path(output.fasta)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Check for local FASTA first
        if params.local_fasta and Path(params.local_fasta).exists():
            print(f"Linking local FASTA: {params.local_fasta}")
            # Copy instead of symlink for portability
            shutil.copy(params.local_fasta, output_path)
            print(f"Copied to {output_path}")

        elif params.proteome_id:
            # Download from UniProt
            url = f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(proteome:{params.proteome_id})"
            print(f"Downloading from UniProt: {params.proteome_id}")
            print(f"URL: {url}")

            urllib.request.urlretrieve(url, output_path)
            print(f"Downloaded to {output_path}")

        else:
            raise ValueError(f"No FASTA source configured for organism: {wildcards.organism}")


rule download_contaminants:
    """
    Download cRAP contaminants database.
    """
    output:
        fasta="resources/fasta/contaminants.fasta"
    params:
        url=config.get("contaminants", {}).get("url", "https://www.thegpm.org/crap/crap.fasta")
    log:
        "logs/fasta/contaminants.log"
    shell:
        """
        curl -L -o {output.fasta} {params.url} 2>&1 | tee {log}
        """


def get_contaminants_input(wildcards):
    """Get contaminants input - empty if disabled or organism already has them."""
    if not config.get("contaminants", {}).get("enabled", True):
        return []
    organism = get_organism_for_dataset(wildcards.accession)
    if organism_includes_contaminants(organism):
        return []  # Organism's local FASTA already includes contaminants
    return "resources/fasta/contaminants.fasta"


rule prepare_search_database:
    """
    Prepare complete search database for a dataset.

    Combines organism FASTA with contaminants (unless organism FASTA already includes them).
    """
    input:
        organism_fasta=lambda wildcards: f"resources/fasta/{get_organism_for_dataset(wildcards.accession)}.fasta",
        contaminants=get_contaminants_input
    output:
        fasta="resources/fasta/search_db/{accession}.fasta"
    log:
        "logs/fasta/search_db/{accession}.log"
    run:
        from pathlib import Path

        output_path = Path(output.fasta)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        print(f"Creating search database for {wildcards.accession}")

        with open(output_path, 'w') as out_f:
            # Add organism FASTA
            print(f"  Adding organism FASTA: {input.organism_fasta}")
            with open(input.organism_fasta) as in_f:
                out_f.write(in_f.read())

            # Add contaminants if enabled
            if input.contaminants:
                print(f"  Adding contaminants: {input.contaminants}")
                with open(input.contaminants) as in_f:
                    out_f.write("\n")
                    out_f.write(in_f.read())

        # Count sequences
        n_seqs = sum(1 for line in open(output_path) if line.startswith('>'))
        print(f"  Total sequences: {n_seqs}")


rule list_fasta_databases:
    """
    List available FASTA databases.
    """
    run:
        from pathlib import Path

        fasta_dir = Path("resources/fasta")
        print("\n=== Available FASTA Databases ===")

        if fasta_dir.exists():
            for f in sorted(fasta_dir.glob("*.fasta")):
                n_seqs = sum(1 for line in open(f) if line.startswith('>'))
                size_mb = f.stat().st_size / (1024 * 1024)
                print(f"  {f.name}: {n_seqs:,} sequences ({size_mb:.1f} MB)")

            search_db = fasta_dir / "search_db"
            if search_db.exists():
                print("\n  Search databases:")
                for f in sorted(search_db.glob("*.fasta")):
                    n_seqs = sum(1 for line in open(f) if line.startswith('>'))
                    print(f"    {f.name}: {n_seqs:,} sequences")
        else:
            print("  No FASTA databases found.")

        print("=================================\n")
