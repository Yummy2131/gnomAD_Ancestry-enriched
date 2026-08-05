import os
import argparse

def delete_duplicates(directory, prefix_len=29, dry_run=False):
    seen = {}
    deleted = []

    for fname in sorted(os.listdir(directory)):
        path = os.path.join(directory, fname)

        # skip directories
        if not os.path.isfile(path):
            continue

        prefix = fname[:prefix_len]

        if prefix in seen:
            if dry_run:
                print(f"[DRY RUN] Would delete: {path}")
            else:
                os.remove(path)
                print(f"Deleted: {path}")
            deleted.append(path)
        else:
            seen[prefix] = path

    print(f"\nSummary:")
    print(f"Kept files: {len(seen)}")
    print(f"Deleted files: {len(deleted)}")


def main():
    parser = argparse.ArgumentParser(
        description="Delete duplicate files based on first N characters of filename"
    )
    parser.add_argument(
        "directory",
        help="Directory containing files"
    )
    parser.add_argument(
        "--prefix-len",
        type=int,
        default=29,
        help="Number of leading characters to compare (default: 29)"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be deleted without deleting"
    )

    args = parser.parse_args()

    delete_duplicates(
        directory=args.directory,
        prefix_len=args.prefix_len,
        dry_run=args.dry_run
    )


if __name__ == "__main__":
    main()
