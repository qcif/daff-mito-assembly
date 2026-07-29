---
name: execute
description: Execute a task from the ./tasks directory
---

- Read ./tasks/$1_*.md and carry out the specified task.
- If anything is unclear, stop and ask for clarification in the task.
- If the task requires the use of a tool or library, ensure that a container is
  available for it. Run the container with docker to confirm that the tool
  works as expected, and that inputs, outputs and arguments are as expected.
  This should be light-touch, and could be as simple as running the tool's
  --help command.
- If any deviations or discoveries are made as the task progresses, write them
  to the bottom of the task file in an "Outcomes" section.
- When the task is complete, ensure that any checklists in the task file are
  checked off, and move the task file to ./tasks/completed/.
- If task execution resulted in new tasks being created, write them to
  ./tasks/todo.md.
