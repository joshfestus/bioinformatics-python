# Solutions

One program per chapter, holding the answers to that chapter's
exercises. Each prints what the exercise asks for and the reasoning
behind it, so the answer is worked rather than asserted.

```
cd solutions
uv run ch01.py
```

Run them from this folder. They read the chapter folders beside them, so
a chapter whose data is fetched rather than committed needs
`data/fetch.py` run first, and the program says so when something is
missing.

A few exercises need the network. Those catch the failure and print the
answer anyway, so being offline costs you the live numbers and nothing
else.

Where an exercise asks you to work something out by hand, the program
gives the hand answer first and then checks it with code. That is the
order the book teaches, and it is the order worth keeping.
