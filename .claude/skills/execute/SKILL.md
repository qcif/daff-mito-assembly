---
name: execute
description: Execute a task from the ./tasks directory
---

Read ./tasks/$1_*.md and carry out the specified task.

**Before starting**
- Prompt the user if `git status` shows a lot of uncommitted
  changes - they should probably be committed before starting.
- Consider whether the task has unmet preconditions, or whether it should be
  deferred until another task is complete.

**When executing the task**
- If anything is unclear, stop and ask for clarification in the task.
- If the task requires the use of a tool or library, ensure that a container is
  available for it. Run the container with docker to confirm that the tool
  works as expected, and that inputs, outputs and arguments are as expected.
  This should be light-touch, and could be as simple as running the tool's
  --help command.
- Deviations from the task brief require approval from the user before
  proceeding.
- As the task progresses, write any outcomes to the bottom of the task file in
  an "Outcomes" section, especially any deviations from the task definition, or
  any new information that comes to light.
- When the task is complete, ensure that any checklists in the task file are
  checked off, and move the task file to ./tasks/completed/.
- Ensure that any markdown links to the task file are updated to point to the
  new location in ./tasks/completed/.
- If task execution resulted in new tasks being created, write them to
  ./tasks/todo.md.
- If modifying spec, do not create markdown links to task files, just say
  "task nn_task_name.md". This prevents broken links when tasks are moved or
  renamed.
