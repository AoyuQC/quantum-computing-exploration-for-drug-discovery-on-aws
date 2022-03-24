当您在亚马逊云科技基础设施上构建解决方案时，安全责任由您和亚马逊云科技共同承担。此[责任共担模型](https://aws.amazon.com/compliance/shared-responsibility-model/)减少了您的操作负担，这是由于亚马逊云科技操作、管理和控制组件，包括主机操作系统、虚拟化层以及服务运行所在设施的物理安全性。有关亚马逊云科技安全的更多信息，请访问亚马逊云科技[云安全](http://aws.amazon.com/security/)。

## 安全最佳实践

药物发现量子计算解决方案在设计时考虑了安全最佳实践。但是，解决方案的安全性会根据您的具体用例而有所不同，有时添加额外的安全措施会增加解决方案的成本。以下是增强此解决方案在生产环境中的安全性的建议。

### IAM角色

亚马逊云科技身份和访问管理（IAM）角色允许客户为亚马逊云科技上的服务和用户分配细粒度访问策略和权限。此解决方案创建IAM角色，这些角色授予解决方案各组件间的访问权限。

### 安全组

此解决方案中创建的安全组旨在控制和隔离各组件间的网络流量。我们建议您检查安全组，并在部署启动并运行后根据需要进一步限制访问。

### Amazon Braket安全设置

由于在解决方案中使用了Amazon Braket，请参考与其相关的[安全措施说明](https://docs.aws.amazon.com/braket/latest/developerguide/security.html)。