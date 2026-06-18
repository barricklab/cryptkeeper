class FakeLogger:
    """An empty logger that does nothing to prevent cryptkeeper from crashing"""

    def __init__(self):
        pass

    def info(self, message):
        pass

    def debug(self, message):
        pass

    def warning(self, message):
        pass

    def error(self, message):
        pass

    def critical(self, message):
        pass


class delay_iterator:
    def __init__(self, iterable):
        self.iterable = iterable
        self.delayed_items = []

    def delay(self, item):
        self.delayed_items.append(item)

    def is_delayed(self):
        if self.iterable:
            return False
        elif self.delayed_items:
            return True
        else:
            return False

    def __iter__(self):
        return self

    def __next__(self):
        if self.iterable:
            return self.iterable.pop(0), False
        elif self.delayed_items:
            return self.delayed_items.pop(0), True
        else:
            raise StopIteration

    def __len__(self):
        return len(self.iterable) + len(self.delayed_items)

    def __repr__(self):
        return f"delay_iterator(iterable={self.iterable}, delayed_items={self.delayed_items})"

    def __str__(self):
        return f"delay_iterator(iterable={self.iterable}, delayed_items={self.delayed_items})"

    def __bool__(self):
        return bool(self.iterable or self.delayed_items)

    def __getitem__(self, key):
        return (
            self.iterable[key]
            if key < len(self.iterable)
            else self.delayed_items[key - len(self.iterable)]
        )


# EOF
