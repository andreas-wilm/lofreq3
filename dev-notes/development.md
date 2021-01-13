https://nvie.com/posts/a-successful-git-branching-model/
he central repo holds two main branches with an infinite lifetime:

master
develop

When the source code in the develop branch reaches a stable point and is ready to be released, all of the changes should be merged back into master somehow and then tagged with a release number.

Creating a feature branch 
When starting work on a new feature, branch off from the develop branch.

$ git checkout -b myfeature develop


## Testing

- tests should mainly be done in module, except complex stuff
- for non-exported functions use `testblock` in `when isMainModule`, which is automatically tested in tests/moduletest.nim
- FIXME log running tests that need data
