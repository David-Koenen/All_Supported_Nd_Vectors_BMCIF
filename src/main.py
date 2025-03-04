import argparse
import os
from run_tests import run_instance_from_file, create_and_run_instance, run_best_tests, run_worst_tests


def main():
    parser = argparse.ArgumentParser(description="Run bi-objective minimum cost flow problem instances.")
    parser.add_argument("--method", type=str, required=True, choices=[
        "run_instance_from_file", "create_and_run_instance", "run_best_tests", "run_worst_tests"],
                        help="Method to execute.")
    parser.add_argument("--nodes", type=int, help="Number of nodes in the graph.")
    parser.add_argument("--arcs", type=int, help="Number of arcs in the graph.")
    parser.add_argument("--number", type=int, help="Instance number (for run_instance_from_file).")
    parser.add_argument("--seed", type=int, help="Random seed (for create_and_run_instance).")
    parser.add_argument("--m1", action='store_true', help="Activate first method.")
    parser.add_argument("--m2", action='store_true', help="Activate second method.")
    parser.add_argument("--m3", action='store_true', help="Activate third method.")
    parser.add_argument("--m4", action='store_true', help="Activate fourth method.")
    parser.add_argument("--output", type=str, default="output/results.txt", help="Output file path.")
    
    args = parser.parse_args()
    
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    if args.method == "run_instance_from_file":
        if None in (args.nodes, args.arcs, args.number):
            raise ValueError("--nodes, --arcs, and --number are required for run_instance_from_file")
        run_instance_from_file(args.nodes, args.arcs, args.number, args.m1, args.m2, args.m3, args.m4)
    
    elif args.method == "create_and_run_instance":
        if None in (args.nodes, args.arcs, args.seed):
            raise ValueError("--nodes, --arcs, and --seed are required for create_and_run_instance")
        create_and_run_instance(args.nodes, args.arcs, args.seed, args.m1, args.m2, args.m3, args.m4)
    
    elif args.method == "run_best_tests":
        run_best_tests()
    
    elif args.method == "run_worst_tests":
        run_worst_tests()
    
    print(f"Results saved to {args.output} and also added to output\Tests_Results.csv")


if __name__ == "__main__":
    main()
