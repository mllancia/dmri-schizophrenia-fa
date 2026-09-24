# Git routine (PyCharm)

Keep this file in each project, or somewhere you can find it.

---

## The mental model

Your work lives in three places. Two commands move it along:

```
your folder  --[ commit ]-->  local history  --[ push ]-->  GitHub
```

- **Commit** = save a snapshot with a note about what changed. Local only.
- **Push** = upload those snapshots to GitHub.

Commit often, push when you want it backed up or shared.

---

## Daily loop — two shortcuts

| Do this | Shortcut | What it means |
|---|---|---|
| Commit | **⌘K** | Snapshot your changes |
| Push | **⌘⇧K** | Send to GitHub |
| Pull | **⌘T** | Get changes from GitHub (only if working on two machines) |
| See history | **⌘9** | The Git tool window |

In the ⌘K dialog: tick the files to include, write a message, click
**Commit and Push** to do both at once.

**Write messages in the present tense, saying what changed:**

- good: `Add motion QC with pre-specified thresholds`
- good: `Exclude sub-10855 for orientation artifact`
- bad: `update`, `stuff`, `asdf`

Future-you reads these when trying to find when something broke.

---

## THE ONE RULE

**Create `.gitignore` before connecting a project to GitHub.**

This is the only step that can cause real trouble. A large file committed by
accident is painful to remove afterwards, and GitHub rejects anything over
100 MB — so the push fails after a long wait.

Right-click the project root → New → File → `.gitignore`, and paste:

```
.venv/
venv/
__pycache__/
*.pyc
.idea/
.DS_Store

data/
raw/
*.nii.gz
*.csv.gz
*.zip
```

Adjust the bottom section to whatever holds data in that project. The rule of
thumb: **if you can re-download it or regenerate it, ignore it.** Code,
configuration, small results tables and write-ups go in. Data and installed
packages stay out.

---

## Starting a new project (once per project)

1. Create `.gitignore` first (see above).
2. **Git → GitHub → Share Project on GitHub**
3. Name it, choose Private or Public, click Share.
4. PyCharm shows the files it will add. **Read that list.** Everything should
   be code, text, or a small results file. If you see anything large or any
   data file, cancel and fix `.gitignore`.
5. Done. From here it's just ⌘K and ⌘⇧K.

If the GitHub menu is missing: **PyCharm → Settings → Version Control →
GitHub → +** to sign in.

If Git isn't enabled at all: **VCS → Enable Version Control Integration →
Git**.

---

## Checking what you're about to send

Before a first push, or any push you're unsure about:

- **⌘9** opens the Git window. The **Log** tab shows your commits; the
  **Commit** tab shows what's staged.
- Files listed in green are new, blue are modified.
- Nothing large or data-shaped should appear. That's the whole check.

---

## When something goes wrong

**"I committed something I shouldn't have" (and haven't pushed)**
⌘9 → Log → right-click the commit → **Undo Commit**. Your files stay; the
snapshot goes. Fix `.gitignore`, commit again.

**"Push was rejected"**
Someone or something changed GitHub since you last pulled. Press **⌘T** to
pull, resolve anything PyCharm flags, then push again.

**"It's asking for a password"**
GitHub stopped accepting passwords. Use a personal access token: GitHub →
Settings → Developer settings → Personal access tokens → generate one with
`repo` scope, and paste it where the password goes. Easier: sign in through
PyCharm's GitHub settings once and it handles this.

**"I want to check GitHub actually has it"**
Open the repo in a browser. If your latest commit message is at the top, the
push worked. Don't trust the local state alone.

---

## Terminal equivalents

Only needed if PyCharm isn't available. Same three ideas.

```
git status              what has changed
git add -A              stage everything not ignored
git commit -m "message" snapshot it
git push                send to GitHub
git pull                fetch from GitHub
git log --oneline       recent commits
```

---

## What belongs in a repository

**Yes:** scripts, notebooks, `README.md`, `requirements.txt`, `.gitignore`,
small results tables, figures, write-ups, decision logs.

**No:** raw data, derived image volumes, virtual environments, anything over
a few MB, credentials or API keys.

The test: if someone cloned this repo and ran your scripts, could they
reproduce the work? If yes, it has everything it needs. Data comes from
wherever you originally got it, and your `README.md` should say where.
