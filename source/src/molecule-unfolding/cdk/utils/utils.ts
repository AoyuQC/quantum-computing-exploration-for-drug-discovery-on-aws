import * as cdk from 'aws-cdk-lib'

import {
    aws_iam as iam
} from 'aws-cdk-lib'

import {
    aws_kms as kms
} from 'aws-cdk-lib'

import {
    aws_ec2 as ec2
} from 'aws-cdk-lib'

import {
    IAspect,
} from 'aws-cdk-lib';

import {
    IConstruct
} from 'constructs'


export class ChangePublicSubnet implements IAspect {
    visit(node: IConstruct): void {
        if (node instanceof ec2.CfnSubnet && node.mapPublicIpOnLaunch) {
            node.addPropertyOverride('MapPublicIpOnLaunch', false)
        }
    }
}

export class AddCfnNag implements IAspect {
    visit(node: IConstruct): void {
        if (

            node.node.path.endsWith('/Custom::S3AutoDeleteObjectsCustomResourceProvider/Handler') ||
            node.node.path.endsWith('/BucketNotificationsHandler050a0587b7544547bf325f094a3db834/Resource')
        ) {
            (node as cdk.CfnResource).addMetadata('cfn_nag', {
                rules_to_suppress: [{
                        id: 'W58',
                        reason: 'the lambda is auto generated by CDK',
                    },
                    {
                        id: 'W89',
                        reason: 'the lambda is auto generated by CDK',
                    },
                    {
                        id: 'W92',
                        reason: 'the lambda is auto generated by CDK',
                    }
                ],
            });
        } else if (
            node.node.path.endsWith('/AggResultLambda/Resource') ||
            node.node.path.endsWith('/TaskParametersLambda/Resource') ||
            node.node.path.endsWith('/DeviceAvailableCheckLambda/Resource') ||
            node.node.path.endsWith('/WaitForTokenLambda/Resource') ||
            node.node.path.endsWith('/BraketTaskEventHandler/ParseBraketResultLambda/Resource')
        ) {
            (node as cdk.CfnResource).addMetadata('cfn_nag', {
                rules_to_suppress: [{
                    id: 'W58',
                    reason: 'the lambda already have the cloudwatch permission',
                }, ],
            });
        } else if (node.node.path.endsWith('/hpcBatchJobRole/DefaultPolicy/Resource') ||
            node.node.path.endsWith('/qcBatchJobRole/DefaultPolicy/Resource') ||
            node.node.path.endsWith('/createModelBatchJobRole/DefaultPolicy/Resource') ||
            node.node.path.endsWith('/batchExecutionRole/DefaultPolicy/Resource') ||
            node.node.path.endsWith('/TaskParametersLambdaRole/DefaultPolicy/Resource') ||
            node.node.path.endsWith('/DeviceAvailableCheckLambdaRole/DefaultPolicy/Resource') ||
            node.node.path.endsWith('/ParseBraketResultLambdaRole/DefaultPolicy/Resource') ||
            node.node.path.endsWith('/AggResultLambdaRole/DefaultPolicy/Resource') ||
            node.node.path.endsWith('/WaitForTokenLambdaRole/DefaultPolicy/Resource') ||
            node.node.path.endsWith('/MolUnfNotebook/NotebookRole/DefaultPolicy/Resource') ||
            node.node.path.endsWith('/BucketNotificationsHandler050a0587b7544547bf325f094a3db834/Role/DefaultPolicy/Resource')
        ) {
            (node as cdk.CfnResource).addMetadata('cfn_nag', {
                rules_to_suppress: [{
                    id: 'W12',
                    reason: 'some permissions are not resource-level permissions',
                }, ],
            });
        } else if (node.node.path.endsWith('/HPCStateMachine/Role/DefaultPolicy/Resource') ||
            node.node.path.endsWith('/QCStateMachine/Role/DefaultPolicy/Resource')
        ) {
            (node as cdk.CfnResource).addMetadata('cfn_nag', {
                rules_to_suppress: [{
                    id: 'W12',
                    reason: 'the policy about log group is generated by CDK',
                }, ],
            });
        } else if (
            node.node.path.endsWith('/BenchmarkStateMachine/Role/DefaultPolicy/Resource') ||
            node.node.path.endsWith('/RunHPCAndQCStateMachine/Role/DefaultPolicy/Resource') ||
            node.node.path.endsWith('/QCDeviceStateMachine/Role/DefaultPolicy/Resource')
        ) {
            (node as cdk.CfnResource).addMetadata('cfn_nag', {
                rules_to_suppress: [{
                        id: 'W12',
                        reason: 'the policy about log group is generated by CDK',
                    },
                    {
                        id: 'W76',
                        reason: 'The policy is generated automatically by CDK',
                    },
                ],
            });
        } else if (node.node.path.endsWith('/AccessLogS3Bucket/Resource')) {
            (node as cdk.CfnResource).addMetadata('cfn_nag', {
                rules_to_suppress: [{
                    id: 'W35',
                    reason: 'this is access log bucket ',
                }, ],
            });
        } else if (node.node.path.endsWith('/batchSg/Resource') ||
            node.node.path.endsWith('/lambdaSg/Resource')
        ) {
            (node as cdk.CfnResource).addMetadata('cfn_nag', {
                rules_to_suppress: [{
                    id: 'W5',
                    reason: 'cidr open to world on egress',
                }],
            });
        }
    }
}

export function grantKmsKeyPerm(key: kms.IKey, logGroupName ? : string): void {
    key.addToResourcePolicy(new iam.PolicyStatement({
        principals: [new iam.ServicePrincipal('logs.amazonaws.com')],
        actions: [
            'kms:Encrypt*',
            'kms:ReEncrypt*',
            'kms:Decrypt*',
            'kms:GenerateDataKey*',
            'kms:Describe*',
        ],
        resources: [
            '*',
        ],
        conditions: {
            ArnLike: {
                'kms:EncryptionContext:aws:logs:arn': cdk.Arn.format({
                    service: 'logs',
                    resource: 'log-group',
                    resourceName: logGroupName ? logGroupName : '*',
                    arnFormat: cdk.ArnFormat.COLON_RESOURCE_NAME,
                }, cdk.Stack.of(key)),
            },
        },
    }));
}