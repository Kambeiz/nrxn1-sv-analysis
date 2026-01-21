#!/usr/bin/env python3
"""AWS Batch job submission utilities for NRXN1 analysis."""

import boto3
import json
from typing import Dict, Any, Optional
from pathlib import Path
import time


class AWSBatchRunner:
    """Submit and manage AWS Batch jobs for NRXN1 analysis."""
    
    def __init__(
        self,
        job_queue: str = "nrxn1-analysis-queue",
        job_definition: str = "nrxn1-pipeline",
        region: str = "us-east-1"
    ):
        """Initialize AWS Batch runner.
        
        Args:
            job_queue: Name of the AWS Batch job queue.
            job_definition: Name of the job definition.
            region: AWS region.
        """
        self.job_queue = job_queue
        self.job_definition = job_definition
        self.region = region
        self.batch_client = boto3.client("batch", region_name=region)
        self.s3_client = boto3.client("s3", region_name=region)
    
    def submit_job(
        self,
        input_s3_path: str,
        output_s3_path: str,
        job_name: Optional[str] = None,
        vcpus: int = 4,
        memory: int = 16000,
        parameters: Optional[Dict[str, str]] = None
    ) -> str:
        """Submit a job to AWS Batch.
        
        Args:
            input_s3_path: S3 path to input VCF file.
            output_s3_path: S3 path for output results.
            job_name: Optional job name.
            vcpus: Number of vCPUs.
            memory: Memory in MB.
            parameters: Additional job parameters.
        
        Returns:
            Job ID.
        """
        if job_name is None:
            job_name = f"nrxn1-analysis-{int(time.time())}"
        
        container_overrides = {
            "vcpus": vcpus,
            "memory": memory,
            "environment": [
                {"name": "INPUT_S3_PATH", "value": input_s3_path},
                {"name": "OUTPUT_S3_PATH", "value": output_s3_path}
            ],
            "command": [
                "python", "main.py",
                "--input", "/data/input.vcf",
                "--output", "/data/output",
                "--mode", "full"
            ]
        }
        
        response = self.batch_client.submit_job(
            jobName=job_name,
            jobQueue=self.job_queue,
            jobDefinition=self.job_definition,
            containerOverrides=container_overrides,
            parameters=parameters or {}
        )
        
        return response["jobId"]
    
    def get_job_status(self, job_id: str) -> Dict[str, Any]:
        """Get status of a submitted job.
        
        Args:
            job_id: Job ID to check.
        
        Returns:
            Job status information.
        """
        response = self.batch_client.describe_jobs(jobs=[job_id])
        
        if response["jobs"]:
            job = response["jobs"][0]
            return {
                "job_id": job["jobId"],
                "job_name": job["jobName"],
                "status": job["status"],
                "created_at": job.get("createdAt"),
                "started_at": job.get("startedAt"),
                "stopped_at": job.get("stoppedAt"),
                "status_reason": job.get("statusReason", "")
            }
        
        return {"job_id": job_id, "status": "NOT_FOUND"}
    
    def wait_for_completion(
        self,
        job_id: str,
        poll_interval: int = 30,
        timeout: int = 3600
    ) -> str:
        """Wait for job completion.
        
        Args:
            job_id: Job ID to wait for.
            poll_interval: Seconds between status checks.
            timeout: Maximum wait time in seconds.
        
        Returns:
            Final job status.
        """
        start_time = time.time()
        
        while True:
            status = self.get_job_status(job_id)
            current_status = status["status"]
            
            if current_status in ["SUCCEEDED", "FAILED"]:
                return current_status
            
            if time.time() - start_time > timeout:
                return "TIMEOUT"
            
            time.sleep(poll_interval)
    
    def download_results(
        self,
        s3_output_path: str,
        local_output_dir: Path
    ) -> None:
        """Download results from S3.
        
        Args:
            s3_output_path: S3 path to results.
            local_output_dir: Local directory to download to.
        """
        bucket, prefix = self._parse_s3_path(s3_output_path)
        
        local_output_dir.mkdir(parents=True, exist_ok=True)
        
        paginator = self.s3_client.get_paginator("list_objects_v2")
        
        for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
            for obj in page.get("Contents", []):
                key = obj["Key"]
                local_path = local_output_dir / Path(key).relative_to(prefix)
                local_path.parent.mkdir(parents=True, exist_ok=True)
                
                self.s3_client.download_file(bucket, key, str(local_path))
    
    def _parse_s3_path(self, s3_path: str) -> tuple:
        """Parse S3 path into bucket and key."""
        path = s3_path.replace("s3://", "")
        parts = path.split("/", 1)
        bucket = parts[0]
        key = parts[1] if len(parts) > 1 else ""
        return bucket, key


def create_job_definition() -> Dict[str, Any]:
    """Create AWS Batch job definition template."""
    return {
        "jobDefinitionName": "nrxn1-pipeline",
        "type": "container",
        "containerProperties": {
            "image": "your-ecr-repo/nrxn1-analysis:latest",
            "vcpus": 4,
            "memory": 16000,
            "command": ["python", "main.py"],
            "jobRoleArn": "arn:aws:iam::ACCOUNT_ID:role/BatchJobRole",
            "volumes": [
                {
                    "name": "data",
                    "host": {"sourcePath": "/data"}
                }
            ],
            "mountPoints": [
                {
                    "sourceVolume": "data",
                    "containerPath": "/data",
                    "readOnly": False
                }
            ],
            "environment": [
                {"name": "AWS_DEFAULT_REGION", "value": "us-east-1"}
            ]
        },
        "retryStrategy": {"attempts": 2},
        "timeout": {"attemptDurationSeconds": 43200}
    }


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Submit NRXN1 analysis to AWS Batch")
    parser.add_argument("--input", required=True, help="S3 path to input VCF")
    parser.add_argument("--output", required=True, help="S3 path for output")
    parser.add_argument("--wait", action="store_true", help="Wait for completion")
    
    args = parser.parse_args()
    
    runner = AWSBatchRunner()
    job_id = runner.submit_job(args.input, args.output)
    print(f"Submitted job: {job_id}")
    
    if args.wait:
        print("Waiting for completion...")
        status = runner.wait_for_completion(job_id)
        print(f"Final status: {status}")
