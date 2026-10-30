# 1.用户及文件信息配置

## 配置用户名和邮箱等信息

```C
git config --global user.name "Kris Wu"
git config --global user.email "你的GitHub邮箱"
```

查看用户名和邮箱

```
git config --global --list
```

## 配置推送地址

```C
git remote set-url origin https://github.com/holyfuxkingshxt/git_test.git
```

清理错误代理

```C
git config --global --unset http.proxy
git config --global --unset https.proxy
```

设置代理

```C
git config --global http.proxy http://127.0.0.1:7890
git config --global https.proxy http://127.0.0.1:7890
//proxy为端口号
```

恢复或者设置地址

```
git remote set-url origin https://github.com/holyfuxkingshxt/git_test.git
```

# 2.文件上传

## 2.1 git add 暂存文件（告诉 Git 哪些文件要上传）

```C
# 只上传单个文件
git add README.md

# 上传当前目录所有修改、新增文件（最常用）
git add .
```

### 对应影响

1. 把你修改 / 新建的文件存入**暂存区**；
2. `git status` 里文件从红色（未暂存）变成绿色（待提交）；
3. **仅存在本地电脑，GitHub 远程仓库完全无变化**；
4. 撤销方式：`git reset HEAD 文件名` 取消暂存。

## 2.2 git commit 本地提交（生成本地版本记录）

```C
git commit -m "填写本次修改说明，必填"
```

### 对应影响

1. 将暂存区的改动打包，生成一条**本地唯一提交记录**（一串哈希编号，如 75240b6）；
2. 生成本地版本快照，就算删文件也能通过这条记录找回；
3. 改动**依旧只存在你电脑本地**，远程仓库看不到任何更新；
4. 可以多次 commit 攒多个版本，最后一次性上传远程。

## 2.3 git push 推送到 GitHub 远程仓库（真正上传到网页）

### 场景 1：分支第一次上传远程（自动创建远程分支）

```C
git push -u origin dev
```

`-u`：绑定本地 dev 分支 ↔ 远程 origin/dev 分支；

远程没有 dev 分支时，自动在 GitHub 新建 dev 分支；

### 场景 2：分支已绑定过上游，后续重复上传（简写）

```C
git push
```

### push 执行后的影响

1. 把你本地所有未同步的 commit 提交，完整上传到 GitHub 远程仓库；
2. 刷新 GitHub 网页，就能看到你修改、新增的文件；
3. 远程分支代码与本地同步；
4. 若远程有别人新提交代码，直接 push 会报错，需要先`git pull`拉取合并；
5. 多人协作时，所有人拉取代码后都能看到你的修改。

# 3.分支管理

## 3.1查看所有分支

```C
# 只看本地分支，* 代表你现在在哪条分支
git branch

# 本地+远程所有分支全部列出
git branch -a
```

```matlab
#1. 查看当前分支全部文件
git ls-files
#执行后会打印仓库内所有文件名称，都是 main 分支保存的文件。
#2. 查看完整目录树结构（带文件夹层级）
git ls-tree -r main
#-r：递归，显示子文件夹里的文件
#main：指定查看 main 分支，就算切到别的分支也能查 main 的文件
# 3. 只看顶层文件夹 / 文件（不深入子目录）
git ls-tree main
```

# 3.2 创建分支

```C
## 创建并切到新分支
# 新版Git推荐
git switch -c dev

# 兼容旧版Git
git checkout -b dev
#dev为分支名
```

```C
##创建并不切到新分支
git branch dev
```

```
## 在已有分支之间切换
# 切到 dev 分支
git switch dev
# 或者旧写法
git checkout dev
```

## 3.3 示例从创建到上传

```matlab
# 1从main创建并切换到dev分支
git switch -c dev
# 2修改文件，暂存、提交
git add .
git commit -m "新增测试功能"
# 3 第一次推送到 GitHub，推送命令自动生成远程分支
git push -u origin dev
# 切回 main 主分支
git switch main
```

注：git push -u origin dev 	

```C
git push：把本地代码上传到远程仓库
origin：GitHub 仓库的远程别名：{GitHub 仓库名：git_test 本地对这个仓库的代称：origin}
dev：你本地的 dev 分支
    
-u = --set-upstream，设置上游关联；绑定本地 dev ↔ 远程 origin/dev
绑定前：推送必须写完整 git push origin dev，拉取也要 git pull origin dev
绑定后：以后直接简写 git push / git pull，Git 自动识别你当前在 dev，只会同步远程 dev，不用再加分支名
```

## 3.4  删除远程分支

```matlab
git push origin --delete dev
```

# 4 远程拉取

## 4.1 轻量拉取：git fetch（只同步信息，不合并）

```C
$ git fetch origin
```

`fetch` = **拉取远程仓库的最新信息**，不改动本地任何代码

1. 同步远程仓库所有分支清单（比如你网页上新建的 dev、别人新建分支，本地会更新 `remotes/origin/*` 列表）；
2. 下载远程每个分支最新的提交记录，更新本地缓存 `remotes/origin/main`、`remotes/origin/dev`；

## 4.2 完整拉取：git pull（下载 + 自动合并，最常用）

在 main 分支：

```C
git pull
等价于：拉取远程 origin/main 并合并到本地 main
```

当前在 dev 分支：

```C
git pull
```

等价于：拉取远程 origin/dev 并合并到本地 dev

1. 先自动执行一次 `fetch` 下载远程全部更新

2. 自动把远程对应分支的新代码合并到你现在本地分支

# 5 合并

把两条分支的代码整合到一起，解决两边各自修改的差异，最终统一成一份代码。

1. 本地两个分支互相合并（dev 开发代码合并到 main 主干）

2. `git pull` 拉取远程时自动合并远程代码到本地分支

## 场景 1：本地分支合并（dev 合并到 main，最常用）

需求：dev 写完功能，把 dev 代码合并到稳定 main 分支

```
步骤 1：切到要接收代码的目标分支（main）
git switch main
步骤 2：保证本地 main 是最新的（先拉远程，避免冲突）
git pull
步骤 3：执行合并，把 dev 的代码合进当前 main
git merge dev
```

合并后两种结果

### ① 无冲突，自动合并成功

终端提示 `Automatic merge went well`

- main 分支现在同时拥有 main + dev 所有代码；
- 本地合并完成，只是电脑本地生效，还没上传 GitHub；
- 最后执行 `git push` 同步远程 main。

### ② 有代码冲突（远程 / 本地改了同一个文件同一行）

Git 分不清保留哪边代码，合并中断：

1. 终端提示 `Automatic merge failed; fix conflicts and then commit the result.`
2. 打开冲突文件，会看到冲突标记：

```C
<<<<<<< HEAD   # 当前main分支代码
原有文字
=======       # 分割线
dev分支修改的文字
>>>>>>> dev    # 来源分支
```

3. **手动编辑文件**：删掉标记符号，保留你想要的最终内容，保存文件

4. 全部冲突改完后，把文件加入暂存区

```C
git add .
```

5. 完成合并提交

```C
git commit -m "解决冲突，合并dev到main"
```

6. 推送远程

```C
git push
```

放弃本次合并（合并一半不想合了）

冲突处理到一半想取消合并，恢复合并前状态：

```C
git merge --abort
```

## 场景 2：拉取远程代码自动合并（git pull = fetch + merge）

当远程 origin/main 有别人新提交，你本地 main 落后，执行：

```C
git pull
```

1. 先自动 `git fetch` 下载远程最新代码

2. 自动执行 `git merge origin/main` 把远程代码合并到本地 main同样会出现「自动合并成功 / 代码冲突」两种情况，冲突解决方式同上。

## 场景 3：远程 dev 合并到本地 dev（多人协作）

你在 dev 分支开发，同事上传了新代码到远程 dev：

```C
git switch dev
git pull
```

自动合并远程 origin/dev 到你本地 dev。

### tips

合并不会删除原分支：dev 合并到 main 后，dev 分支依旧存在，可继续开发；

想删除合并完不用的分支：

本地删除：`git branch -d dev`

远程删除：`git push origin --delete dev`

推荐工作流：日常在 dev 写代码，稳定后合并到 main，main 永远保留可运行稳定版本。