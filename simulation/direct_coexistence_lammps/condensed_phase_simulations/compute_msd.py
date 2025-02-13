import subprocess

def extract_dimensions(atom_filename):
    """get box dimensions (Lx, Lz) from the specified .atom file"""
    with open(atom_filename, 'r') as file:
        lines = file.readlines()
    Lx = float(lines[5].split()[1]) - float(lines[5].split()[0])  # extract Lx value
    Lz = float(lines[7].split()[1]) - float(lines[7].split()[0])  # extract Lz value
    return Lx, Lz

def extract_timestep_size(seqid):
    """get the timestep size from the temp-cycle.in file for a given seqID"""
    filepath = f"{seqid}/temp-cycle.in"
    with open(filepath, 'r') as file:
        lines = file.readlines()
    timestep_size = float(lines[121].split()[1])
    return timestep_size

def process_msd_files(seqid, trial_num, timestep_size):
    """process MSD files for this seq and trial # """
    msd_file = f"msd_temp-cycle-trial{trial_num}.dat"
    atom_file = f"temp-cycle-trial{trial_num}-12.atom"

    Lx, Lz = extract_dimensions(atom_file)
    subprocess.run([
        'python3', '../get_starting_phase_keyfiles/msd_plot.py', 
        seqid, msd_file, str(Lz), str(Lx), str(timestep_size)
    ])

def main():
    """Run MSD analysis for all sequences listed in seqs-ps.lst"""
    trial_numbers = [1]  # currently using one trial
    with open('seqs-ps.lst', 'r') as file:
        for line in file:
            seqid = line.split()[0]
            timestep_size = extract_timestep_size(seqid)
            for trial_num in trial_numbers:
                process_msd_files(seqid, trial_num, timestep_size)

if __name__ == "__main__":
    main()
