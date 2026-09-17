import io
import os
import polars as pl


def load_outlet_post(
        file_path: str
    ) -> pl.DataFrame:
    """ Load and parse outlet.post file using polars.

    Parameters
    ----------
    file_path : str
        The absolute or relative path to the outlet.post file.

    Returns
    -------
    pl.DataFrame
        The parsed data as a Polars DataFrame.
    """
    with open(file_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    if not lines:
        return pl.DataFrame()

    header_line = lines[0].strip()

    if header_line.startswith("#"):
        header_line = header_line[1:].strip()

    header_clean = header_line.replace("(", " ").replace(")", " ")
    header_cols = header_clean.split()

    data_rows = []

    for line in lines[1:]:
        l_strip = line.strip()
        if not l_strip or l_strip.startswith("#"):
            continue

        cleaned = l_strip.replace("(", " ").replace(")", " ")
        data_rows.append("\t".join(cleaned.split()))

    if not data_rows:
        return pl.DataFrame()

    first_row_cols = data_rows[0].split("\t")
    num_cols = len(first_row_cols)

    if "Y1..YN" in header_cols:
        idx = header_cols.index("Y1..YN")
        num_ys = num_cols - len(header_cols) + 1
        new_ys = [f"Y{i+1}" for i in range(num_ys)]
        header_cols = header_cols[:idx] + new_ys + header_cols[idx+1:]

    elif len(header_cols) < num_cols:
        for i in range(len(header_cols), num_cols):
            header_cols.append(f"col_{i}")

    elif len(header_cols) > num_cols:
        header_cols = header_cols[:num_cols]

    csv_content = "\t".join(header_cols) + "\n" + "\n".join(data_rows)

    return pl.read_csv(io.StringIO(csv_content), separator="\t")


if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    outlet_post_path = os.path.join(
        script_dir,
        "postProcessing",
        "lagrangian",
        "cloud",
        "funcPatchPostProcessing",
        "1000",
        "outlet.post",
    )

    if os.path.exists(outlet_post_path):
        df = load_outlet_post(outlet_post_path)
        print(f"Loaded DataFrame with shape: {df.shape}")
        print(df.head(5).to_dicts())

    else:
        print(f"File not found: {outlet_post_path}")
