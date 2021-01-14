https://nvie.com/posts/a-successful-git-branching-model/
the central repo holds two main branches with an infinite lifetime:

- master
- develop

When the source code in the develop branch reaches a stable point and is ready
to be released, all of the changes should be merged back into master somehow
and then tagged with a release number.
o master should always pass all tests

## Creating a feature branch 

When starting work on a new feature, branch off from the develop branch.

$ git checkout -b myfeature develop

When done merge into devel and run all tests (nimble tests and others)
before merging into master (see above)

## Git


    # once feature is complete:
    # on devel
    git merge --squash feature-branch
    git commit
    # run tests

    # later on master
    # run tests
    git merge devel
    git tag -a v3.0.0-rc1 -m "version 3.0.0 release candidate 1"
    git push origin --tags



## Testing

- tests should mainly be done in module, except complex stuff
- for non-exported functions use `testblock` in `when isMainModule`, which is automatically tested in tests/moduletest.nim
- FIXME log running tests that need data
