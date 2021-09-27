import csv
import math

from .config import Config
from .region import Region
from .types import PrimerDirection


class DeepScheme:
    """A scheme with n pools"""

    def __init__(self, fasta: tuple[str, str], cfg: Config) -> None:
        self.fasta = fasta
        self.cfg = cfg

        self.ref = fasta[1]
        self.ref_id = fasta[0]
        max_insert_percent = cfg.amplicon_size_max / cfg.insert_size_max
        self.pools: list[list[Region]] = [
            []
            for _ in range(math.ceil((cfg.coverage / cfg.packing) * max_insert_percent))
        ]
        self.region_size = int(1 / cfg.packing * cfg.amplicon_size_max)
        self.pool_offset = int(self.region_size / len(self.pools))

    def run(self) -> None:
        for n, pool in enumerate(self.pools):
            offset = n * self.pool_offset
            print(f"Pool {n}, offset {offset}")
            for i, start in enumerate(range(offset, len(self.ref), self.region_size)):
                print(f"Region {i} starting at {start}")
                ref_slice = self.ref[
                    start : min((start + self.region_size), len(self.ref))
                ]
                if len(ref_slice) >= self.cfg.amplicon_size_max * 2:
                    region = Region(
                        ref_slice,
                        start,
                        self.cfg,
                    )
                    print(region)
                    pool.append(region)
                else:
                    print("Pool end")

    def __str__(self) -> str:
        return (
            f"DeepScheme: {self.fasta[0]}, region size: {self.region_size} ({self.cfg})"
        )

    def write_primer_bed(self) -> None:
        """Write primer BED file."""
        filepath = self.cfg.output / f"{self.cfg.prefix}.primer.bed"

        rows = []
        ref_id = self.ref_id
        primer_name_template = "POOL_{}_REGION_{}_{}"

        for pool_num, pool in enumerate(self.pools):
            for region_num, region in enumerate(pool):
                for direction in PrimerDirection:
                    name = primer_name_template.format(
                        pool_num, region_num, direction.name
                    )
                    p = getattr(region.top_pair, direction.name.lower())
                    start = p.start if p.direction == PrimerDirection.FORWARD else p.end
                    end = p.end if p.direction == PrimerDirection.FORWARD else p.start
                    rows.append(
                        [
                            ref_id,
                            start,
                            end + 1,  # BED is half-open
                            name,
                            pool_num,
                            direction.value,
                            p.seq,
                        ]
                    )

            with filepath.open("w") as fh:
                cw = csv.writer(fh, delimiter="\t")
                cw.writerows(rows)

    def write_insert_bed(self) -> None:
        """Write insert BED file."""
        filepath = self.cfg.output / f"{self.cfg.prefix}.insert.bed"

        rows = []
        ref_id = self.ref_id
        insert_name_template = "POOL_{}_REGION_{}"

        for pool_num, pool in enumerate(self.pools):
            for region_num, region in enumerate(pool):
                if region.insert:
                    name = insert_name_template.format(pool_num, region_num)
                    insert = region.insert
                    rows.append(
                        [
                            ref_id,
                            insert.start,
                            insert.end + 1,  # BED is half-open
                            name,
                            "+",
                        ]
                    )

            with filepath.open("w") as fh:
                cw = csv.writer(fh, delimiter="\t")
                cw.writerows(rows)
