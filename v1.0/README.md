
# TripinRNA

## Introduction

This project contains scripts to identify MALAT1 or NEAT1-like triplex structures in entire transcriptomes. The code is written in two versions: 
- `no_gap.py`: For finding triplex structures without gaps in the Hoogsteen strand.
- `with_gap.py`: For finding triplex structures with gaps in the Hoogsteen strand.

The scripts are flexible and can be used on any type of transcriptome from any species. They have been tested on human protein-coding and lncRNA transcriptomes. The output is a CSV file containing the identified triplex structures.

**Note**: The raw files in the `raw_backup` directory were used to generate the output for our manuscript. The new versions (v1.0) are more organized and potentially faster.

## Prerequisites

To run these scripts, you need:

- Python 3.6 or later
- The following Python packages:
  - pandas
  - gzip
  - shutil
  - re
  - typing
  - concurrent
  - datetime

## Installation

1. **Clone the repository**:
   ```bash
   git clone https://github.com/yourusername/TripinRNA.git
   ```
2. **Install dependencies**:
   - Option 1: Using `pip`:
     ```bash
     pip install -r requirements.txt
     ```
   - Option 2: Using `conda`:
     ```bash
     conda create --name triplex_env python=3.8
     conda activate triplex_env
     pip install -r requirements.txt
     ```

## Usage

### Running on Local System

1. **Without Gap**:
   ```bash
   python src/no_gap.py
   ```

2. **With Gap**:
   ```bash
   python src/with_gap.py
   ```

### Running on HPC with Slurm

1. **Without Gap**:
   ```bash
   sbatch scripts/no_gap_slurm.sh
   ```

2. **With Gap**:
   ```bash
   sbatch scripts/with_gap_slurm.sh
   ```

## Output

The scripts output a CSV file containing the identified triplex structures. The format is as follows:

| Column Name              | Description                                     |
|--------------------------|-------------------------------------------------|
| Hoogsteen Strand         | The Hoogsteen strand sequence                    |
| H Start Index            | Start index of the Hoogsteen strand             |
| H End Index              | End index of the Hoogsteen strand               |
| Crick Strand             | The Crick strand sequence                        |
| C Start Index            | Start index of the Crick strand                 |
| C End Index              | End index of the Crick strand                   |
| Watson Strand            | The Watson strand sequence                       |
| W Start Index            | Start index of the Watson strand                |
| W End Index              | End index of the Watson strand                  |
| Upper Loop Sequence      | The upper loop sequence                          |
| Upper Loop Start Index   | Start index of the upper loop                   |
| Upper Loop End Index     | End index of the upper loop                     |
| Lower Loop Sequence      | The lower loop sequence                          |
| Lower Loop Start Index   | Start index of the lower loop                   |
| Lower Loop End Index     | End index of the lower loop                     |
| Lower Stem (downstream) Sequence | The downstream lower stem sequence     |
| Lower Stem (downstream) Start Index | Start index of the downstream lower stem |
| Lower Stem (downstream) End Index | End index of the downstream lower stem  |
| Lower Stem (upstream) Sequence | The upstream lower stem sequence         |
| Lower Stem (upstream) Start Index | Start index of the upstream lower stem |
| Lower Stem (upstream) End Index | End index of the upstream lower stem     |
| Pseudoknot Sequence      | The pseudoknot sequence                         |
| Pseudoknot Start Index   | Start index of the pseudoknot                   |
| Pseudoknot End Index     | End index of the pseudoknot                     |
| Pseudoknot Status        | Indicates if a pseudoknot is present            |

## License

This project is licensed under the MIT License. See the LICENSE file for details.
