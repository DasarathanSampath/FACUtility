#!/usr/bin/env python3
"""
Process ParaView CSV files containing points on circular elbow slices.

For each CSV file, the script:
1. Fits the actual slice plane from Points:0, Points:1, Points:2.
2. Fits a circle in that plane.
3. Calculates circumferential angle using:
       0°   = longest-radius side (extrados)
       90°  = positive side selected by POSITIVE_90_REFERENCE
       180° = shortest-radius side (intrados)
       270° = opposite side
4. Sorts points by angle.
5. Plots TLPatch against circumferential angle.
6. Interpolates TLPatch at 0°, 45°, ..., 315°.
7. Calculates TLPatch and yPlus statistics in 45° sectors.
8. Writes processed CSV files and combined summary tables.

Required packages:
    pip install numpy pandas matplotlib
"""

from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# USER SETTINGS
# ---------------------------------------------------------------------------

# Folder containing the ParaView CSV files.
INPUT_FOLDER = Path("comparisonData/simulated/case3_slices")

# File pattern. Examples:
# "*.csv"       -> process every CSV file
# "slice_*.csv" -> process files such as slice_4.csv, slice_30.csv, etc.
FILE_PATTERN = "*.csv"

# Output folder.
OUTPUT_FOLDER = Path("processed_slices")

# Coordinate columns exported by ParaView.
POINT_COLUMNS = ["Points:0", "Points:1", "Points:2"]

# Field to plot and tabulate.
VALUE_COLUMN = "TLPatch"

# Additional field for binned statistics.
YPLUS_COLUMN = "yPlus"

# Coordinates of the elbow/bend centre.
#
# This point defines the shortest-radius direction:
# fitted slice centre -> BEND_CENTER = 180° direction.
#
# Replace these values with the centre of curvature of your elbow.
BEND_CENTER = np.array([0.0, 0.0, 0.0], dtype=float)

# Reference direction used to choose the +90° side.
# For an elbow lying in the global X-Y plane, global +Z is normally suitable.
# Change to [0, 0, -1] if 90° and 270° appear reversed.
POSITIVE_90_REFERENCE = np.array([0.0, 0.0, 1.0], dtype=float)

# Requested circumferential positions.
TARGET_ANGLES = np.arange(0.0, 360.0, 45.0)

# True: interpolate between adjacent circumferential points.
# False: select the nearest available point.
USE_INTERPOLATION = True

# Plot all slices in one figure as well as individual figures.
CREATE_COMBINED_PLOT = True

NORMALIZED_VALUE_COLUMN = "TLPatch_Normalized"

# ---------------------------------------------------------------------------
# GEOMETRY FUNCTIONS
# ---------------------------------------------------------------------------

def find_overall_maximum(
    csv_files: list[Path],
    value_column: str,
) -> float:
    """
    Find the maximum TLPatch value across all input CSV files.
    """
    overall_maximum = -np.inf

    for csv_path in csv_files:
        dataframe = pd.read_csv(
            csv_path,
            usecols=[value_column],
        )

        values = pd.to_numeric(
            dataframe[value_column],
            errors="coerce",
        ).to_numpy(dtype=float)

        valid_values = values[np.isfinite(values)]

        if len(valid_values) == 0:
            continue

        file_maximum = np.max(valid_values)
        overall_maximum = max(overall_maximum, file_maximum)

    if not np.isfinite(overall_maximum):
        raise ValueError(
            f"No valid values were found in column '{value_column}'."
        )

    if np.isclose(overall_maximum, 0.0):
        raise ValueError(
            "The overall maximum is zero. Normalization cannot be performed."
        )

    return float(overall_maximum)

def calculate_binned_statistics(
    angles_deg: np.ndarray,
    values: np.ndarray,
    target_angles: np.ndarray,
    half_width: float = 22.5,
    value_name: str = "TLPatch",
) -> pd.DataFrame:
    """
    Calculate count, minimum, maximum, and mean within angular bins.

    Each bin is centered on a target angle and extends:
        target_angle ± half_width

    Circular distance is used, so the bin centered at 0 degrees includes:
        337.5° <= angle < 360°
        0° <= angle <= 22.5°

    Parameters
    ----------
    angles_deg
        Circumferential angles in degrees, within 0 to 360.
    values
        TLPatch values corresponding to the angles.
    target_angles
        Centers of the angular bins.
    half_width
        Half-width of each bin in degrees. Default is 22.5°.

    Returns
    -------
    pandas.DataFrame
        Statistics for each target angle.
    """
    angles_deg = np.asarray(angles_deg, dtype=float)
    values = np.asarray(values, dtype=float)

    rows = []

    for target in target_angles:

        # Signed circular difference in the range [-180, 180).
        circular_difference = (
            angles_deg - target + 180.0
        ) % 360.0 - 180.0

        # Use a half-open interval to avoid assigning boundary points
        # to two adjacent bins:
        #     -22.5 <= difference < 22.5
        mask = (
            (circular_difference >= -half_width)
            & (circular_difference < half_width)
            & np.isfinite(values)
        )

        selected_values = values[mask]

        if len(selected_values) == 0:
            minimum = np.nan
            maximum = np.nan
            mean = np.nan
            standard_deviation = np.nan
        else:
            minimum = np.min(selected_values)
            maximum = np.max(selected_values)
            mean = np.mean(selected_values)
            standard_deviation = np.std(
                selected_values,
                ddof=0,
            )

        rows.append(
            {
                "Target_Angle_deg": target,
                "Bin_Start_deg": (target - half_width) % 360.0,
                "Bin_End_deg": (target + half_width) % 360.0,
                "Point_Count": len(selected_values),
                f"{value_name}_Min": minimum,
                f"{value_name}_Max": maximum,
                f"{value_name}_Mean": mean,
                f"{value_name}_Standard_Deviation": standard_deviation,
            }
        )

    return pd.DataFrame(rows)

def normalize(vector: np.ndarray, tolerance: float = 1.0e-12) -> np.ndarray:
    """Return a unit vector."""
    magnitude = np.linalg.norm(vector)
    if magnitude < tolerance:
        raise ValueError("A required direction vector has near-zero magnitude.")
    return vector / magnitude


def fit_plane(points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Fit a plane using singular-value decomposition.

    Returns
    -------
    reference : ndarray, shape (3,)
        Mean point on the plane.
    normal : ndarray, shape (3,)
        Unit plane normal.
    """
    reference = points.mean(axis=0)
    _, _, vh = np.linalg.svd(points - reference, full_matrices=False)
    normal = normalize(vh[-1])
    return reference, normal


def make_slice_basis(
    normal: np.ndarray,
    positive_90_reference: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Create two orthogonal directions in the fitted plane.

    The first direction is only an intermediate fitting basis.
    The final angular zero direction is determined from BEND_CENTER.
    """
    candidate = positive_90_reference - normal * np.dot(
        positive_90_reference, normal
    )

    # Fallback when the supplied reference is almost normal to the slice.
    if np.linalg.norm(candidate) < 1.0e-10:
        fallback_vectors = (
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            np.array([0.0, 0.0, 1.0]),
        )
        for fallback in fallback_vectors:
            candidate = fallback - normal * np.dot(fallback, normal)
            if np.linalg.norm(candidate) >= 1.0e-10:
                break

    basis_v = normalize(candidate)
    basis_u = normalize(np.cross(basis_v, normal))
    return basis_u, basis_v


def fit_circle_2d(u: np.ndarray, v: np.ndarray) -> tuple[float, float, float]:
    """
    Algebraic least-squares circle fit:
        (u - uc)^2 + (v - vc)^2 = radius^2
    """
    matrix = np.column_stack((2.0 * u, 2.0 * v, np.ones(len(u))))
    rhs = u**2 + v**2
    uc, vc, constant = np.linalg.lstsq(matrix, rhs, rcond=None)[0]

    radius_squared = constant + uc**2 + vc**2
    if radius_squared <= 0.0:
        raise ValueError("Circle fit produced a non-positive radius.")

    return float(uc), float(vc), float(np.sqrt(radius_squared))


def calculate_circumferential_angles(
    points: np.ndarray,
    bend_center: np.ndarray,
    positive_90_reference: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float, np.ndarray, float]:
    """
    Calculate circumferential angles.

    Convention
    ----------
    0°   : extrados, opposite the bend centre
    90°  : side selected by POSITIVE_90_REFERENCE
    180° : intrados, toward the bend centre
    270° : opposite side

    Returns
    -------
    angles_deg
    radii
    fitted_radius
    circle_center_3d
    maximum_plane_residual
    """
    reference, normal = fit_plane(points)
    basis_u, basis_v = make_slice_basis(normal, positive_90_reference)

    relative = points - reference
    u = relative @ basis_u
    v = relative @ basis_v

    uc, vc, fitted_radius = fit_circle_2d(u, v)
    circle_center_3d = reference + uc * basis_u + vc * basis_v

    radial_vectors = points - circle_center_3d
    radii = np.linalg.norm(radial_vectors, axis=1)

    # Intrados direction: from section centre toward bend centre.
    inward = bend_center - circle_center_3d
    inward -= normal * np.dot(inward, normal)
    inward = normalize(inward)

    # 0° is the extrados, opposite the inward direction.
    zero_direction = -inward

    # Project the requested +90° reference into the plane and make it
    # perpendicular to the zero direction.
    ninety_direction = positive_90_reference.copy()
    ninety_direction -= normal * np.dot(ninety_direction, normal)
    ninety_direction -= zero_direction * np.dot(
        ninety_direction, zero_direction
    )
    ninety_direction = normalize(ninety_direction)

    # Make the 90° direction follow the supplied reference rather than its
    # opposite.
    projected_reference = (
        positive_90_reference
        - normal * np.dot(positive_90_reference, normal)
    )
    if np.dot(ninety_direction, projected_reference) < 0.0:
        ninety_direction *= -1.0

    x_component = radial_vectors @ zero_direction
    y_component = radial_vectors @ ninety_direction

    angles_deg = np.mod(
        np.degrees(np.arctan2(y_component, x_component)),
        360.0,
    )

    plane_residual = np.abs((points - reference) @ normal)
    maximum_plane_residual = float(plane_residual.max())

    return (
        angles_deg,
        radii,
        fitted_radius,
        circle_center_3d,
        maximum_plane_residual,
    )


# ---------------------------------------------------------------------------
# DATA EXTRACTION FUNCTIONS
# ---------------------------------------------------------------------------

def periodic_interpolate(
    angles_deg: np.ndarray,
    values: np.ndarray,
    target_angles: np.ndarray,
) -> np.ndarray:
    """
    One-dimensional periodic linear interpolation over 0° to 360°.

    Duplicate angles are averaged before interpolation.
    """
    clean = pd.DataFrame(
        {"angle": angles_deg, "value": values}
    ).dropna()

    # Average duplicate or numerically identical angular positions.
    clean["angle_key"] = clean["angle"].round(10)
    grouped = (
        clean.groupby("angle_key", as_index=False)
        .agg(angle=("angle", "mean"), value=("value", "mean"))
        .sort_values("angle")
    )

    x = grouped["angle"].to_numpy(dtype=float)
    y = grouped["value"].to_numpy(dtype=float)

    if len(x) < 2:
        raise ValueError("At least two valid angular points are required.")

    # Extend at both ends to enforce periodic interpolation.
    x_periodic = np.concatenate(([x[-1] - 360.0], x, [x[0] + 360.0]))
    y_periodic = np.concatenate(([y[-1]], y, [y[0]]))

    return np.interp(target_angles, x_periodic, y_periodic)


def nearest_values(
    angles_deg: np.ndarray,
    values: np.ndarray,
    target_angles: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return values and actual angles of nearest circumferential points."""
    selected_values = []
    selected_angles = []

    for target in target_angles:
        circular_difference = np.abs(
            (angles_deg - target + 180.0) % 360.0 - 180.0
        )
        index = int(np.nanargmin(circular_difference))
        selected_values.append(values[index])
        selected_angles.append(angles_deg[index])

    return np.asarray(selected_values), np.asarray(selected_angles)


def natural_sort_key(path: Path):
    """Sort filenames containing numbers in human order."""
    return [
        int(token) if token.isdigit() else token.lower()
        for token in re.split(r"(\d+)", path.stem)
    ]


# ---------------------------------------------------------------------------
# FILE PROCESSING
# ---------------------------------------------------------------------------

def process_csv(
    csv_path: Path,
    overall_maximum: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    """Process one ParaView CSV file."""
    dataframe = pd.read_csv(csv_path)

    dataframe[VALUE_COLUMN] = pd.to_numeric(
        dataframe[VALUE_COLUMN],
        errors="coerce",
    )

    dataframe[YPLUS_COLUMN] = pd.to_numeric(
        dataframe[YPLUS_COLUMN],
        errors="coerce",
    )

    dataframe[NORMALIZED_VALUE_COLUMN] = (
        dataframe[VALUE_COLUMN] / overall_maximum
    )

    original_values = dataframe[
        VALUE_COLUMN
    ].to_numpy(dtype=float)

    normalized_values = dataframe[
        NORMALIZED_VALUE_COLUMN
    ].to_numpy(dtype=float)

    yplus_values = dataframe[
        YPLUS_COLUMN
    ].to_numpy(dtype=float)

    missing = [
        column
        for column in [*POINT_COLUMNS, VALUE_COLUMN, YPLUS_COLUMN]
        if column not in dataframe.columns
    ]
    if missing:
        raise KeyError(
            f"{csv_path.name}: missing required columns: {missing}"
        )

    points = dataframe[POINT_COLUMNS].to_numpy(dtype=float)

    (
        angles,
        radii,
        fitted_radius,
        fitted_center,
        max_plane_residual,
    ) = calculate_circumferential_angles(
        points=points,
        bend_center=BEND_CENTER,
        positive_90_reference=POSITIVE_90_REFERENCE,
    )

    dataframe["Circumferential_Angle_deg"] = angles
    dataframe["Radius_from_fitted_center"] = radii
    dataframe = dataframe.sort_values(
        "Circumferential_Angle_deg"
    ).reset_index(drop=True)

    values = dataframe[
        NORMALIZED_VALUE_COLUMN
    ].to_numpy(dtype=float)
    sorted_angles = dataframe["Circumferential_Angle_deg"].to_numpy(
        dtype=float
    )

    statistics_table = calculate_binned_statistics(
        angles_deg=sorted_angles,
        values=values,
        target_angles=TARGET_ANGLES,
        half_width=22.5,
    )

    original_statistics = calculate_binned_statistics(
        angles_deg=sorted_angles,
        values=original_values,
        target_angles=TARGET_ANGLES,
        half_width=22.5,
        value_name="TLPatch",
    )

    normalized_statistics = calculate_binned_statistics(
        angles_deg=sorted_angles,
        values=normalized_values,
        target_angles=TARGET_ANGLES,
        half_width=22.5,
        value_name="TLPatch_Normalized",
    )

    yplus_statistics = calculate_binned_statistics(
        angles_deg=sorted_angles,
        values=yplus_values,
        target_angles=TARGET_ANGLES,
        half_width=22.5,
        value_name="yPlus",
    )

    statistics_table = original_statistics.merge(
        normalized_statistics[
            [
                "Target_Angle_deg",
                "TLPatch_Normalized_Min",
                "TLPatch_Normalized_Max",
                "TLPatch_Normalized_Mean",
                "TLPatch_Normalized_Standard_Deviation",
            ]
        ],
        on="Target_Angle_deg",
        how="left",
    )

    statistics_table = statistics_table.merge(
        yplus_statistics[
            [
                "Target_Angle_deg",
                "yPlus_Min",
                "yPlus_Max",
                "yPlus_Mean",
                "yPlus_Standard_Deviation",
            ]
        ],
        on="Target_Angle_deg",
        how="left",
    )

    statistics_table.insert(
        0,
        "Slice",
        csv_path.stem,
    )

    if USE_INTERPOLATION:
        selected_values = periodic_interpolate(
            sorted_angles,
            values,
            TARGET_ANGLES,
        )
        selected_actual_angles = TARGET_ANGLES.copy()
        selection_method = "periodic linear interpolation"
    else:
        selected_values, selected_actual_angles = nearest_values(
            sorted_angles,
            values,
            TARGET_ANGLES,
        )
        selection_method = "nearest point"

    table = pd.DataFrame(
        {
            "Slice": csv_path.stem,
            "Target_Angle_deg": TARGET_ANGLES,
            "Actual_or_Interpolated_Angle_deg": selected_actual_angles,
            VALUE_COLUMN: selected_values,
        }
    )

    metadata = {
        "Slice": csv_path.stem,
        "Number_of_points": len(dataframe),
        "Circle_center_X": fitted_center[0],
        "Circle_center_Y": fitted_center[1],
        "Circle_center_Z": fitted_center[2],
        "Fitted_radius": fitted_radius,
        "Mean_point_radius": radii.mean(),
        "Maximum_radius_deviation": np.max(
            np.abs(radii - fitted_radius)
        ),
        "Maximum_plane_residual": max_plane_residual,
        "Selection_method": selection_method,
    }

    return dataframe, table, statistics_table, metadata


def save_individual_plot(
    processed: pd.DataFrame,
    slice_name: str,
    output_path: Path,
) -> None:
    """Create a line plot for one slice."""
    figure, axis = plt.subplots(figsize=(8.5, 5.0))

    angle_values = processed[
    "Circumferential_Angle_deg"
    ].to_numpy(dtype=float)

    tlpatch_values = processed[
        NORMALIZED_VALUE_COLUMN
    ].to_numpy(dtype=float)

    axis.plot(
        angle_values,
        tlpatch_values,
        linewidth=1.5,
    )

    axis.scatter(
        angle_values,
        tlpatch_values,
        s=12,
    )

    axis.set_xlim(0.0, 360.0)
    axis.set_xticks(np.arange(0.0, 361.0, 45.0))
    axis.set_xlabel("Circumferential angle (°)")
    axis.set_ylabel(
        r"Normalized $T_L$ ($T_L$ / overall maximum)"
    )
    axis.set_ylim(0.0, 1.05)
    axis.set_title(f"{VALUE_COLUMN} around circumference — {slice_name}")
    axis.grid(True, alpha=0.3)
    figure.tight_layout()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    """Process all matching slice files."""
    OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)

    csv_files = sorted(
        INPUT_FOLDER.glob(FILE_PATTERN),
        key=natural_sort_key,
    )

    if not csv_files:
        raise FileNotFoundError(
            f"No files matching '{FILE_PATTERN}' were found in "
            f"'{INPUT_FOLDER.resolve()}'."
        )
    
    overall_tlpatch_maximum = find_overall_maximum(
        csv_files=csv_files,
        value_column=VALUE_COLUMN,
    )

    print(
        f"Overall maximum {VALUE_COLUMN}: "
        f"{overall_tlpatch_maximum:.10g}"
    )

    all_tables = []
    metadata_rows = []
    combined_series = []
    all_statistics_tables = []

    for csv_path in csv_files:
        print(f"Processing: {csv_path.name}")

        (
            processed,
            target_table,
            statistics_table,
            metadata,
        ) = process_csv(
            csv_path=csv_path,
            overall_maximum=overall_tlpatch_maximum,
        )

        processed_csv = (
            OUTPUT_FOLDER / f"{csv_path.stem}_with_angles.csv"
        )
        table_csv = (
            OUTPUT_FOLDER / f"{csv_path.stem}_{VALUE_COLUMN}_at_45deg.csv"
        )
        plot_file = (
            OUTPUT_FOLDER / f"{csv_path.stem}_{VALUE_COLUMN}_plot.png"
        )

        processed.to_csv(processed_csv, index=False)
        target_table.to_csv(table_csv, index=False)
        save_individual_plot(processed, csv_path.stem, plot_file)

        all_tables.append(target_table)
        metadata_rows.append(metadata)
        combined_series.append(
            (
                csv_path.stem,
                processed["Circumferential_Angle_deg"].to_numpy(),
                processed[NORMALIZED_VALUE_COLUMN].to_numpy(),
            )
        )
        statistics_csv = (
            OUTPUT_FOLDER
            / f"{csv_path.stem}_{VALUE_COLUMN}_binned_statistics.csv"
        )

        statistics_table.to_csv(
            statistics_csv,
            index=False,
        )

        all_statistics_tables.append(statistics_table)

    combined_statistics = pd.concat(
        all_statistics_tables,
        ignore_index=True,
    )

    combined_statistics.to_csv(
        OUTPUT_FOLDER
        / "all_planes_TLPatch_binned_statistics_long.csv",
        index=False,
    )

    combined_table = pd.concat(all_tables, ignore_index=True)
    combined_table.to_csv(
        OUTPUT_FOLDER / f"all_slices_{VALUE_COLUMN}_at_45deg_long.csv",
        index=False,
    )

    # ---------------------------------------------------------------
    # Create table:
    # rows    = planes
    # columns = circumferential angles
    # values  = TLPatch
    # ---------------------------------------------------------------

    wide_table = combined_table.pivot(
        index="Slice",
        columns="Target_Angle_deg",
        values=VALUE_COLUMN,
    )

    # Ensure the columns are always in the required order.
    wide_table = wide_table.reindex(columns=TARGET_ANGLES)

    # Rename slice filenames into readable plane names.
    def format_plane_name(slice_name):
        name = str(slice_name)

        # Examples:
        # slice_A -> Plane A
        # plane_A -> Plane A
        # Plane_A -> Plane A
        name = name.replace("slice_", "")
        name = name.replace("Slice_", "")
        name = name.replace("plane_", "")
        name = name.replace("Plane_", "")
        name = name.replace("_", " ")

        return f"Plane {name.strip()}"


    wide_table.index = [
        format_plane_name(name)
        for name in wide_table.index
    ]

    # Rename angle columns.
    wide_table.columns = [
        f"{int(angle)}°"
        for angle in wide_table.columns
    ]

    # Give the row-index column a suitable name.
    wide_table.index.name = "Plane"

    # Save as CSV.
    wide_table.to_csv(
        OUTPUT_FOLDER / "TLPatch_by_plane_and_angle.csv"
    )

    print("\nTLPatch table:")
    print(wide_table)

    pd.DataFrame(metadata_rows).to_csv(
        OUTPUT_FOLDER / "slice_geometry_check.csv",
        index=False,
    )

    if CREATE_COMBINED_PLOT:
        figure, axis = plt.subplots(figsize=(9.0, 5.5))

        for slice_name, angles, values in combined_series:
            axis.plot(
                angles,
                values,
                linewidth=1.8,
                label=slice_name,
            )

        TITLE_FONT_SIZE = 18
        AXIS_LABEL_FONT_SIZE = 16
        TICK_FONT_SIZE = 14
        LEGEND_FONT_SIZE = 13

        axis.set_xlim(0.0, 360.0)
        axis.set_xticks(np.arange(0.0, 361.0, 45.0))
        axis.set_xlabel("Circumferential angle (°)", fontsize=AXIS_LABEL_FONT_SIZE)
        axis.set_ylabel("$T_{L}$, $T_{L}$/overall maximum", fontsize=AXIS_LABEL_FONT_SIZE)
        axis.set_title(
            "Predicted Thickness Loss ($T_{L}$)",
            fontsize=TITLE_FONT_SIZE
        )
        axis.tick_params(
            axis="both",
            which="major",
            labelsize=TICK_FONT_SIZE,
        )
        axis.grid(True, alpha=0.3)
        axis.legend(fontsize=LEGEND_FONT_SIZE)
        figure.tight_layout()
        figure.savefig(
            OUTPUT_FOLDER / f"all_slices_{VALUE_COLUMN}_plot.png",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close(figure)

    print(f"\nCompleted. Results are in: {OUTPUT_FOLDER.resolve()}")

    def format_plane_name(slice_name):
        name = str(slice_name)

        name = name.replace("slice_", "")
        name = name.replace("Slice_", "")
        name = name.replace("plane_", "")
        name = name.replace("Plane_", "")
        name = name.replace("_", " ")

        return f"Plane {name.strip()}"

    for statistic_column, output_name in [
        ("TLPatch_Min", "TLPatch_min_by_plane_and_angle.csv"),
        ("TLPatch_Max", "TLPatch_max_by_plane_and_angle.csv"),
        ("TLPatch_Mean", "TLPatch_mean_by_plane_and_angle.csv"),
        ("yPlus_Min", "yPlus_min_by_plane_and_angle.csv"),
        ("yPlus_Max", "yPlus_max_by_plane_and_angle.csv"),
        ("yPlus_Mean", "yPlus_mean_by_plane_and_angle.csv"),
        ("Point_Count", "point_count_by_plane_and_angle.csv"),
    ]:

        summary_table = combined_statistics.pivot(
            index="Slice",
            columns="Target_Angle_deg",
            values=statistic_column,
        )

        summary_table = summary_table.reindex(
            columns=TARGET_ANGLES
        )

        summary_table.index = [
            format_plane_name(name)
            for name in summary_table.index
        ]

        summary_table.columns = [
            f"{int(angle)}°"
            for angle in summary_table.columns
        ]

        summary_table.index.name = "Plane"

        summary_table.to_csv(
            OUTPUT_FOLDER / output_name
        )

    combined_wide_statistics = combined_statistics.pivot(
        index="Slice",
        columns="Target_Angle_deg",
        values=[
            "TLPatch_Min",
            "TLPatch_Max",
            "TLPatch_Mean",
        ],
    )

    # Put angle first and statistic second.
    combined_wide_statistics = combined_wide_statistics.swaplevel(
        0,
        1,
        axis=1,
    )

    combined_wide_statistics = combined_wide_statistics.reindex(
        columns=pd.MultiIndex.from_product(
            [
                TARGET_ANGLES,
                [
                    "TLPatch_Min",
                    "TLPatch_Max",
                    "TLPatch_Mean",
                ],
            ]
        )
    )

    combined_wide_statistics.index = [
        format_plane_name(name)
        for name in combined_wide_statistics.index
    ]

    combined_wide_statistics.index.name = "Plane"

    combined_wide_statistics.columns = [
        f"{int(angle)}°_{statistic.replace('TLPatch_', '')}"
        for angle, statistic in combined_wide_statistics.columns
    ]

    combined_wide_statistics.to_csv(
        OUTPUT_FOLDER
        / "TLPatch_min_max_mean_by_plane_and_angle.csv"
    )

    # One row per plane with paired TLmean and yPlusMean columns:
    # TLmean_0°, yPlusMean_0°, TLmean_45°, yPlusMean_45°, ...
    paired_mean_table = combined_statistics.pivot(
        index="Slice",
        columns="Target_Angle_deg",
        values=["TLPatch_Mean", "yPlus_Mean"],
    )

    paired_mean_table = paired_mean_table.swaplevel(0, 1, axis=1)
    paired_mean_table = paired_mean_table.reindex(
        columns=pd.MultiIndex.from_product(
            [TARGET_ANGLES, ["TLPatch_Mean", "yPlus_Mean"]]
        )
    )

    paired_mean_table.index = [
        format_plane_name(name)
        for name in paired_mean_table.index
    ]
    paired_mean_table.index.name = "Plane"

    paired_mean_table.columns = [
        (
            f"TLmean_{int(angle)}°"
            if statistic == "TLPatch_Mean"
            else f"yPlusMean_{int(angle)}°"
        )
        for angle, statistic in paired_mean_table.columns
    ]

    paired_mean_table.to_csv(
        OUTPUT_FOLDER / "TLmean_yPlusMean_by_plane_and_angle.csv"
    )
    
if __name__ == "__main__":
    main()
