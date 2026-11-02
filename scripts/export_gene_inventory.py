"""Export complete provider gene/transcript inventories to typed Parquet.

Run with an existing pyarrow environment. Coordinates stay 1-based inclusive;
source attributes and parent relationships become nested fields. SQLite retains
the full exon/CDS/UTR feature inventory. Existing MyGene exports are untouched.
"""
from pathlib import Path
import argparse
import json
import sqlite3
import pyarrow as pa
import pyarrow.parquet as pq
from extension_support import ROOT, save, start, finish, sources
from notebook_support import read_json, file_digest


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inventory',type=Path)
    args=parser.parse_args();config=read_json(ROOT/'config/gene_inventory.json')
    inventory=args.inventory or ROOT/config['output_directory'];output=inventory/'parquet'
    dbpath=inventory/'dog_ros.sqlite';receipt=read_json(inventory/'receipt.private.json')
    if file_digest(dbpath)!=receipt['outputs']['dog_ros.sqlite']:raise ValueError('SQLite inventory changed')
    plan=dict(sources=sources([Path(__file__),ROOT/'scripts/extension_support.py',ROOT/'scripts/notebook_support.py',dbpath,
                              inventory/'receipt.private.json']),pyarrow=pa.__version__)
    if not start(output,plan):return
    db=sqlite3.connect(f'file:{dbpath.as_posix()}?mode=ro',uri=True);db.row_factory=sqlite3.Row
    counts={}
    for table in ['genes','transcripts']:
        records=[]
        for row in db.execute('SELECT * FROM '+table+' ORDER BY record_key'):
            record=dict(row)
            for key in list(record):
                if key.endswith('_json'):
                    record[key[:-5]]=json.loads(record.pop(key))
            records.append(record)
        data=pa.Table.from_pylist(records)
        data=data.replace_schema_metadata({b'assembly':b'ROS_Cfam_1.0',b'coordinate_convention':b'1-based inclusive',
            b'identity':b'Provider records remain distinct; record_key is unique',b'source_sqlite_sha256':file_digest(dbpath).encode()})
        path=output/(table+'.parquet');pq.write_table(data,path,compression='zstd')
        if not pq.read_table(path).equals(data):raise ValueError('Parquet round-trip mismatch')
        counts[table]=len(records);print(table,len(records),'verified Parquet rows',flush=True)
    db.close();finish(output,plan,dict(status='typed_parquet_complete',rows=counts,round_trip_verified=True,
        original_mygene_resource_modified=False,provider_records_merged=False))


if __name__=='__main__':main()
