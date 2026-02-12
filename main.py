import argparse
from pathlib import Path
import sys

src_path = Path(__file__).resolve().parent / "src"
sys.path.append(str(src_path))
from controler import Controller

def main():
    parser = argparse.ArgumentParser( description = 'One Cell Wonder')
    parser.add_argument("folderpath", help="path to the folder with: rules.txt and initial_cell.txt")
    parser.add_argument("-x","--X", help="width",type=int,default=100)
    parser.add_argument("-y","--Y", help="height",type=int,default=100)
    parser.add_argument("--death", help="This tag allow a cell to die.", action='store_true', default=False)

    args=parser.parse_args()

    folder = Path(args.folderpath)
    if not folder.exists() or not folder.is_dir():
        print(f"The folder '{folder}' does not exist or is not a directory.")
        sys.exit(1)

    # Check for required files
    tryFiles = True #usage of non default files
    if tryFiles == True:
        rules_file = folder / "rules.txt"
        initial_file = folder / "initial_cell.txt"
    else:
        rules_file = folder / "rules.txt"
        initial_file = folder / "initial_cell.txt"
    
    missing_files = [f.name for f in (rules_file, initial_file) if not f.exists()]
    if missing_files:
        print(f"Error: The following required file(s) are missing: {', '.join(missing_files)}")
        sys.exit(1)

    c = Controller(args.X,args.Y,initial_file, rules_file)#,args.death)

    
if __name__ == '__main__':
    main()