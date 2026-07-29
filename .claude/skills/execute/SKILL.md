---
name: execute
description: Execute a task from the ./tasks directory
---

- Read ./tasks/$1*.md and carry out the specified task.
- If anything is unclear, stop and ask for clarification in the task.
- If any deviations or discoveries are made as the task progresses, write them to the bottom of the task file in an "Outcomes" section.
- When the task is complete, ensure that any checklists in the task file are checked off, and move the task file to ./tasks/completed/.
